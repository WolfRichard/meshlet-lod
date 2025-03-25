#include <directx/d3dx12.h>
#include <ViewDependentMeshletLoD.h>

#include <psapi.h>

#include "HLSLnames.h"

#include "../shaders/ViewDependentStructures.fxh"

#include <Application.h>
#include <CommandQueue.h>
#include <Helpers.h>
#include <Window.h>
#include <sstream>

#include <imgui.h>
#include "imgui_impl_win32.h"
#include "imgui_impl_dx12.h"
#include "imgui_internal.h"

#include <wrl.h>
using namespace Microsoft::WRL;

#include <d3dcompiler.h>

//#include "DirectXTex.h"

#include <algorithm> 
#if defined(min)
#undef min
#endif

#if defined(max)
#undef max
#endif

using namespace DirectX;


// Clamp a value between a min and max range.
template<typename T>
constexpr const T& clamp(const T& val, const T& min, const T& max)
{
    return val < min ? min : val > max ? max : val;
}

ViewDependentMeshletLoD::ViewDependentMeshletLoD(const std::wstring& name, int width, int height, bool vSync)
    : super(name, width, height, vSync)
    , m_ScissorRect(CD3DX12_RECT(0, 0, LONG_MAX, LONG_MAX))
    , m_Viewport(CD3DX12_VIEWPORT(0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height)))
    , m_FoV(90.0)
    , m_ContentLoaded(false)
{
    m_scene.init(m_model_file_path);
}

void ViewDependentMeshletLoD::UpdateBufferResource(
    ComPtr<ID3D12GraphicsCommandList7> commandList,
    ID3D12Resource** pDestinationResource,
    ID3D12Resource** pIntermediateResource,
    size_t numElements, size_t elementSize, const void* bufferData,
    D3D12_RESOURCE_FLAGS flags)
{
    auto device = Application::Get().GetDevice();

    size_t bufferSize = numElements * elementSize;

    // Create a committed resource for the GPU resource in a default heap.
    ThrowIfFailed(device->CreateCommittedResource(
        &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
        D3D12_HEAP_FLAG_NONE,
        &CD3DX12_RESOURCE_DESC::Buffer(bufferSize, flags),
        D3D12_RESOURCE_STATE_COMMON,
        nullptr,
        IID_PPV_ARGS(pDestinationResource)));

    // Create an committed resource for the upload.
    if (bufferData)
    {
        ThrowIfFailed(device->CreateCommittedResource(
            &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
            D3D12_HEAP_FLAG_NONE,
            &CD3DX12_RESOURCE_DESC::Buffer(bufferSize),
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(pIntermediateResource)));

        D3D12_SUBRESOURCE_DATA subresourceData = {};
        subresourceData.pData = bufferData;
        subresourceData.RowPitch = bufferSize;
        subresourceData.SlicePitch = subresourceData.RowPitch;

        UpdateSubresources(commandList.Get(),
            *pDestinationResource, *pIntermediateResource,
            0, 0, 1, &subresourceData);
    }
}

void ViewDependentMeshletLoD::createPSO()
{
    Application::Get().Flush();
    m_PipelineState.Reset();
    auto device = Application::Get().GetDevice();

    struct PipelineStateStream
    {
        CD3DX12_PIPELINE_STATE_STREAM_ROOT_SIGNATURE pRootSignature;
        CD3DX12_PIPELINE_STATE_STREAM_PRIMITIVE_TOPOLOGY PrimitiveTopologyType;
        CD3DX12_PIPELINE_STATE_STREAM_AS AS;
        CD3DX12_PIPELINE_STATE_STREAM_MS MS;
        CD3DX12_PIPELINE_STATE_STREAM_PS PS;
        CD3DX12_PIPELINE_STATE_STREAM_DEPTH_STENCIL_FORMAT DSVFormat;
        CD3DX12_PIPELINE_STATE_STREAM_RENDER_TARGET_FORMATS RTVFormats;
        CD3DX12_PIPELINE_STATE_STREAM_RASTERIZER RasterizerState;
    } pipelineStateStream;

    D3D12_RT_FORMAT_ARRAY rtvFormats = {};
    rtvFormats.NumRenderTargets = 1;
    rtvFormats.RTFormats[0] = DXGI_FORMAT_R8G8B8A8_UNORM;

    // Define the custom rasterizer state to set culling mode
    CD3DX12_RASTERIZER_DESC rasterizerDesc(D3D12_DEFAULT);
    if (!m_backFaceCulling) rasterizerDesc.CullMode = D3D12_CULL_MODE_NONE;
    if (m_wireframe) rasterizerDesc.FillMode = D3D12_FILL_MODE_WIREFRAME;

    pipelineStateStream.pRootSignature = m_RootSignature.Get();
    pipelineStateStream.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
    pipelineStateStream.MS = CD3DX12_SHADER_BYTECODE(m_meshShaderBlob.Get());
    pipelineStateStream.AS = CD3DX12_SHADER_BYTECODE(m_taskShaderBlob.Get());
    pipelineStateStream.PS = CD3DX12_SHADER_BYTECODE(m_pixelShaderBlob.Get());
    pipelineStateStream.DSVFormat = DXGI_FORMAT_D32_FLOAT;
    pipelineStateStream.RTVFormats = rtvFormats;
    pipelineStateStream.RasterizerState = rasterizerDesc;

    D3D12_PIPELINE_STATE_STREAM_DESC pipelineStateStreamDesc = {
        sizeof(PipelineStateStream), &pipelineStateStream
    };
    ThrowIfFailed(device->CreatePipelineState(&pipelineStateStreamDesc, IID_PPV_ARGS(&m_PipelineState)));
}


void ViewDependentMeshletLoD::setupConstantsUploadBuffer()
{
    auto device = Application::Get().GetDevice();

    // constants buffer
    D3D12_HEAP_PROPERTIES heapProps = {};
    heapProps.Type = D3D12_HEAP_TYPE_UPLOAD;

    D3D12_RESOURCE_DESC constBufferDesc = {};
    constBufferDesc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
    constBufferDesc.Width = (sizeof(S_Constants) + 255) & ~255; // allign to 256 bytes
    constBufferDesc.Height = 1;
    constBufferDesc.DepthOrArraySize = 1;
    constBufferDesc.MipLevels = 1;
    constBufferDesc.Format = DXGI_FORMAT_UNKNOWN;
    constBufferDesc.SampleDesc.Count = 1;
    constBufferDesc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
    constBufferDesc.Flags = D3D12_RESOURCE_FLAG_NONE;

    // Create an upload buffer for constants
    device->CreateCommittedResource(
        &heapProps,
        D3D12_HEAP_FLAG_NONE,
        &constBufferDesc,
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(&m_ConstantsBuffer)
    );
    // Map the buffer (persistent mapping)
    m_ConstantsBuffer->Map(0, nullptr, reinterpret_cast<void**>(&m_mappedConstantData));
}


template <typename T>
void ViewDependentMeshletLoD::setupSrvAndBuffer(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                                                D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle, 
                                                D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                                std::vector<T>& cpuBuffer,
                                                Microsoft::WRL::ComPtr<ID3D12Resource>& gpuBuffer,
                                                std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                                Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                                unsigned int descriptorSize, D3D12_RESOURCE_FLAGS bufferFlags)
{
    auto device = Application::Get().GetDevice();
    srvGpuHandle = nextAvailableGpuSrvHandle;

    copyBuffers.push_back(ComPtr<ID3D12Resource>());
    UpdateBufferResource(commandList, gpuBuffer.GetAddressOf(), copyBuffers.back().GetAddressOf(),
        (uint)cpuBuffer.size(), sizeof(T), cpuBuffer.data(), bufferFlags);

    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
    srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    srvDesc.Buffer.FirstElement = 0;
    srvDesc.Buffer.NumElements = (uint)cpuBuffer.size();
    srvDesc.Buffer.StructureByteStride = sizeof(T);
    srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

    device->CreateShaderResourceView(gpuBuffer.Get(), &srvDesc, nextAvailableCpuSrvHandle);
    nextAvailableCpuSrvHandle.ptr += descriptorSize;
    nextAvailableGpuSrvHandle.ptr += descriptorSize;
}

template <typename T>
void ViewDependentMeshletLoD::setupBindlessSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle, 
                                                         D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle, 
                                                         D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                                         std::vector<std::vector<T>*>& cpuBuffers, 
                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                                         Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                                         unsigned int descriptorSize)
{
    auto device = Application::Get().GetDevice();
    srvGpuHandle = nextAvailableGpuSrvHandle;
    
    for (uint i = 0; i < cpuBuffers.size(); i++)
    {
        std::vector<T>* cpuSingleBuffer = cpuBuffers[i];

        gpuBuffers.push_back(ComPtr<ID3D12Resource>());
        copyBuffers.push_back(ComPtr<ID3D12Resource>());
        
        UpdateBufferResource(commandList, gpuBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)cpuSingleBuffer->size(), sizeof(T), cpuSingleBuffer->data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)cpuSingleBuffer->size();
        srvDesc.Buffer.StructureByteStride = sizeof(T);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(gpuBuffers.back().Get(), &srvDesc, nextAvailableCpuSrvHandle);
        nextAvailableCpuSrvHandle.ptr += descriptorSize;
        nextAvailableGpuSrvHandle.ptr += descriptorSize;
    }
}

// compatible with tight packing of 4 chars into unsigned integers
void ViewDependentMeshletLoD::setupBindlessPrimitiveIndicesSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                                                                         D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle,
                                                                         D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                                                         std::vector<std::vector<unsigned char>*>& primitiveIndicesCpuBuffers,
                                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                                                         Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                                                         unsigned int descriptorSize)
{
    auto device = Application::Get().GetDevice();
    srvGpuHandle = nextAvailableGpuSrvHandle;

    for (uint i = 0; i < primitiveIndicesCpuBuffers.size(); i++)
    {
        std::vector<unsigned char>* cpuSingleBuffer = primitiveIndicesCpuBuffers[i];
        // fill with zeros to allign its size to 32bits
        while (cpuSingleBuffer->size() % 4 != 0)
            cpuSingleBuffer->push_back(0);

        gpuBuffers.push_back(ComPtr<ID3D12Resource>());
        copyBuffers.push_back(ComPtr<ID3D12Resource>());

        UpdateBufferResource(commandList, gpuBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)cpuSingleBuffer->size() / 4, sizeof(uint), cpuSingleBuffer->data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)cpuSingleBuffer->size() / 4;
        srvDesc.Buffer.StructureByteStride = sizeof(uint);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(gpuBuffers.back().Get(), &srvDesc, nextAvailableCpuSrvHandle);
        nextAvailableCpuSrvHandle.ptr += descriptorSize;
        nextAvailableGpuSrvHandle.ptr += descriptorSize;
    }
}



bool ViewDependentMeshletLoD::LoadContent()
{
    auto device = Application::Get().GetDevice();
    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_COPY);
    auto commandList = commandQueue->GetCommandList();

    setupConstantsUploadBuffer();

    D3D12_DESCRIPTOR_HEAP_DESC gpuHeapDesc = {};
    gpuHeapDesc.NumDescriptors = 1                           // single scene objects buffer
                               + 1                           // single work queue buffer
                               + 1                           // single work queue counter buffer
                               + 1                           // single work queue counter clear values buffer
                               + (uint)m_scene.m_mesh_count  // meshlets buffer per unique mesh
                               + (uint)m_scene.m_mesh_count  // meshlet groups buffer per unique mesh
                               + (uint)m_scene.m_mesh_count  // vertex indices buffer per unique mesh
                               + (uint)m_scene.m_mesh_count  // primitive indices buffer per unique mesh
                               + (uint)m_scene.m_mesh_count; // vertex buffer per unique mesh
    gpuHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    gpuHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    gpuHeapDesc.NodeMask = 0;
    device->CreateDescriptorHeap(&gpuHeapDesc, IID_PPV_ARGS(&m_CBV_SRV_UAV_Heap));




    unsigned int descriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
    D3D12_CPU_DESCRIPTOR_HANDLE nextCpuSrvHandle(m_CBV_SRV_UAV_Heap->GetCPUDescriptorHandleForHeapStart());
    D3D12_GPU_DESCRIPTOR_HANDLE nextGpuSrvHandle(m_CBV_SRV_UAV_Heap->GetGPUDescriptorHandleForHeapStart());
    


    std::vector<ComPtr<ID3D12Resource>> copyBuffers;
    
    // vertex-indices buffers
    setupBindlessSrvAndBuffers(m_VertexIndicesSrvHandle, 
                               nextGpuSrvHandle, nextCpuSrvHandle, 
                               m_scene.m_vertex_indices, m_VertexIndicesBuffers, 
                               copyBuffers, commandList, descriptorSize);
    
    // vertex buffers
    setupBindlessSrvAndBuffers(m_VerticesSrvHandle,
                               nextGpuSrvHandle, nextCpuSrvHandle,
                               m_scene.m_vertices, m_VertexBuffers,
                               copyBuffers, commandList, descriptorSize);

    // primitive-indices buffers 
    setupBindlessPrimitiveIndicesSrvAndBuffers(m_PrimitiveIndicesSrvHandle,
                                               nextGpuSrvHandle, nextCpuSrvHandle,
                                               m_scene.m_primitive_indices, m_PrimitiveIndicesBuffers,
                                               copyBuffers, commandList, descriptorSize);

    // meshlets buffers
    setupBindlessSrvAndBuffers(m_MeshletsSrvHandle,
                               nextGpuSrvHandle, nextCpuSrvHandle,
                               m_scene.m_meshlets, m_MeshletBuffers,
                               copyBuffers, commandList, descriptorSize);    

    // meshlet groups buffers
    setupBindlessSrvAndBuffers(m_MeshletGroupsSrvHandle,
                               nextGpuSrvHandle, nextCpuSrvHandle,
                               m_scene.m_meshlet_groups, m_MeshletGroupsBuffers,
                               copyBuffers, commandList, descriptorSize);

    // scene objects buffer
    setupSrvAndBuffer(m_ObjectsSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      m_scene.m_scene_objects, m_ObjectsBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_NONE);

    // work queue buffer
    std::vector<S_WorkQueueEntry> temporaryCpuWorkQueueBufferData(WORK_QUEUE_SIZE, { 0, 0 });
    setupSrvAndBuffer(m_WorkQueueSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      temporaryCpuWorkQueueBufferData, m_WorkQueueBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);

    // work queue counters buffer and clear values
    std::vector<uint> temporaryCpuWorkQueueCountersData(2, 0);
    setupSrvAndBuffer(m_WorkQueueCountersSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      temporaryCpuWorkQueueCountersData, m_WorkQueueCountersBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);
    setupSrvAndBuffer(m_WorkQueueCountersClearValuesSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      temporaryCpuWorkQueueCountersData, m_WorkQueueCountersClearValuesBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_NONE);
    
   

    // Create the descriptor heap for the depth-stencil view.
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&m_DSVHeap)));

    // Load shaders.
    ThrowIfFailed(D3DReadFileToBlob(L"PixelShader.cso", &m_pixelShaderBlob));
    ThrowIfFailed(D3DReadFileToBlob(L"ViewDependentMeshShader.cso", &m_meshShaderBlob));
    ThrowIfFailed(D3DReadFileToBlob(L"ViewDependentTaskShader.cso", &m_taskShaderBlob));

    // Create a root signature.
    D3D12_FEATURE_DATA_ROOT_SIGNATURE featureData = {};
    featureData.HighestVersion = D3D_ROOT_SIGNATURE_VERSION_1_1;
    if (FAILED(device->CheckFeatureSupport(D3D12_FEATURE_ROOT_SIGNATURE, &featureData, sizeof(featureData))))
    {
        featureData.HighestVersion = D3D_ROOT_SIGNATURE_VERSION_1_0;
    }

    // Allow input layout and deny unnecessary access to certain pipeline stages.
    D3D12_ROOT_SIGNATURE_FLAGS rootSignatureFlags =
        D3D12_ROOT_SIGNATURE_FLAG_NONE;


    CD3DX12_ROOT_PARAMETER1 rootParameters[9];
    // constants
    rootParameters[0].InitAsConstantBufferView(0, 0);

    // vertices
    CD3DX12_DESCRIPTOR_RANGE1 srvVerticesRange;
    srvVerticesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 2);
    rootParameters[1].InitAsDescriptorTable(1, &srvVerticesRange, D3D12_SHADER_VISIBILITY_MESH);

    // vertex-indices
    CD3DX12_DESCRIPTOR_RANGE1 srvVertexIndicesRange;
    srvVertexIndicesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 3);
    rootParameters[2].InitAsDescriptorTable(1, &srvVertexIndicesRange, D3D12_SHADER_VISIBILITY_MESH);

    // primitive-indices
    CD3DX12_DESCRIPTOR_RANGE1 srvPrimitiveIndicesRange;
    srvPrimitiveIndicesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 4);
    rootParameters[3].InitAsDescriptorTable(1, &srvPrimitiveIndicesRange, D3D12_SHADER_VISIBILITY_MESH);

    // meshlets
    CD3DX12_DESCRIPTOR_RANGE1 srvMeshletsRange;
    srvMeshletsRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 1);
    rootParameters[4].InitAsDescriptorTable(1, &srvMeshletsRange, D3D12_SHADER_VISIBILITY_ALL);

    // meshlet groups
    CD3DX12_DESCRIPTOR_RANGE1 srvMeshletGroupsRange;
    srvMeshletGroupsRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 5);
    rootParameters[5].InitAsDescriptorTable(1, &srvMeshletGroupsRange, D3D12_SHADER_VISIBILITY_AMPLIFICATION);

    // scene objects
    rootParameters[6].InitAsShaderResourceView(0, 0);

    // work queue objects
    rootParameters[7].InitAsUnorderedAccessView(0, 0);

    // work queue counters objects
    rootParameters[8].InitAsUnorderedAccessView(1, 0);

    
    CD3DX12_VERSIONED_ROOT_SIGNATURE_DESC rootSignatureDescription;
    rootSignatureDescription.Init_1_1(_countof(rootParameters), rootParameters, 0, nullptr, rootSignatureFlags);

    // Serialize the root signature.
    ComPtr<ID3DBlob> rootSignatureBlob;
    ComPtr<ID3DBlob> errorBlob;
    ThrowIfFailed(D3DX12SerializeVersionedRootSignature(&rootSignatureDescription,
        featureData.HighestVersion, &rootSignatureBlob, &errorBlob));

    // Create the root signature.
    ThrowIfFailed(device->CreateRootSignature(0, rootSignatureBlob->GetBufferPointer(),
        rootSignatureBlob->GetBufferSize(), IID_PPV_ARGS(&m_RootSignature)));

    createPSO();

    auto fenceValue = commandQueue->ExecuteCommandList(commandList);
    commandQueue->WaitForFenceValue(fenceValue);

    m_ContentLoaded = true;

    initImGui();

    // Resize/Create the depth buffer.
    ResizeDepthBuffer(GetClientWidth(), GetClientHeight());

    return true;
}

void ViewDependentMeshletLoD::ResizeDepthBuffer(int width, int height)
{
    if (m_ContentLoaded)
    {
        // Flush any GPU commands that might be referencing the depth buffer.
        Application::Get().Flush();

        width = std::max(1, width);
        height = std::max(1, height);

        auto device = Application::Get().GetDevice();

        // Resize screen dependent resources.
        // Create a depth buffer.
        D3D12_CLEAR_VALUE optimizedClearValue = {};
        optimizedClearValue.Format = DXGI_FORMAT_D32_FLOAT;
        optimizedClearValue.DepthStencil = { 1.0f, 0 };

        ThrowIfFailed(device->CreateCommittedResource(
            &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
            D3D12_HEAP_FLAG_NONE,
            &CD3DX12_RESOURCE_DESC::Tex2D(DXGI_FORMAT_D32_FLOAT, width, height,
                1, 0, 1, 0, D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL),
            D3D12_RESOURCE_STATE_DEPTH_WRITE,
            &optimizedClearValue,
            IID_PPV_ARGS(&m_DepthBuffer)
        ));

        // Update the depth-stencil view.
        D3D12_DEPTH_STENCIL_VIEW_DESC dsv = {};
        dsv.Format = DXGI_FORMAT_D32_FLOAT;
        dsv.ViewDimension = D3D12_DSV_DIMENSION_TEXTURE2D;
        dsv.Texture2D.MipSlice = 0;
        dsv.Flags = D3D12_DSV_FLAG_NONE;

        device->CreateDepthStencilView(m_DepthBuffer.Get(), &dsv,
            m_DSVHeap->GetCPUDescriptorHandleForHeapStart());
    }
}

void ViewDependentMeshletLoD::OnResize(ResizeEventArgs& e)
{
    if (e.Width != GetClientWidth() || e.Height != GetClientHeight())
    {
        super::OnResize(e);

        m_Viewport = CD3DX12_VIEWPORT(0.0f, 0.0f,
            static_cast<float>(e.Width), static_cast<float>(e.Height));

        ResizeDepthBuffer(e.Width, e.Height);
    }
}

void ViewDependentMeshletLoD::UnloadContent()
{
    m_ContentLoaded = false;

    m_RootSignature.Reset();
    m_PipelineState.Reset();

    // clear & reset model buffers
    m_CBV_SRV_UAV_Heap.Reset();

    m_MeshletBuffers.clear();
    m_VertexIndicesBuffers.clear();
    m_PrimitiveIndicesBuffers.clear();
    m_VertexBuffers.clear();
    m_MeshletGroupsBuffers.clear();
    m_ObjectsBuffer.Reset();
    m_ConstantsBuffer.Reset();
    m_WorkQueueBuffer.Reset();
    m_WorkQueueCountersBuffer.Reset();

    m_commandSignature.Reset();
}

void ViewDependentMeshletLoD::OnUpdate(UpdateEventArgs& e)
{
    static uint64_t frameCount = 0;
    static double totalTime = 0.0;

   
    m_frameTime = e.ElapsedTime;

    super::OnUpdate(e);

    totalTime += e.ElapsedTime;
    m_totalRunTime += e.ElapsedTime;
    frameCount++;

    if (totalTime > 1.0)
    {
        m_fps = frameCount / totalTime;
        //OutputDebugString(MatrixToString(m_scene.animator.m_FinalBoneMatrices[1]).c_str());
        

        char buffer[512];
        sprintf_s(buffer, "FPS: %f\n", m_fps);
        OutputDebugStringA(buffer);
        //OutputDebugStringA(MatrixToString(m_scene.animator.m_CurrentAnimation->m_RootNode.children[1].children[0].transformation).c_str());
        //OutputDebugStringA(MatrixToString(m_scene.animator.m_CurrentAnimation->FindBone(m_scene.animator.m_CurrentAnimation->m_RootNode.children[1].children[0].name)->m_LocalTransform).c_str());
        //OutputDebugStringA(MatrixToString(DirectX::XMMatrixTranslation(0, 1, 0) * DirectX::XMMatrixTranslation(0, 1, 0)).c_str());

        frameCount = 0;
        totalTime = 0.0;
    }


    // Update the view matrix.
    const XMVECTOR eyePosition = XMVectorSet(0, 0, 0, 1);
    const XMVECTOR focusPoint = XMVectorSet(1, 0, 0, 1);
    const XMVECTOR upDirection = XMVectorSet(0, 1, 0, 0);
    m_ViewMatrix = XMMatrixLookAtLH(eyePosition, focusPoint, upDirection);
    
    // Camera Controll
    if (m_freeCamera)
    {
        float4x4 RotationMatrix = XMMatrixTranspose(XMMatrixRotationY(m_CameraYaw + m_autoRotationOffset) * XMMatrixRotationX(m_CameraRoll));

        if (ImGui::IsKeyDown(ImGuiKey_W)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(0.0f, 0.0f, -1.0f, 0.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_A)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(1.0f, 0.0f, 0.0f, 0.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_S)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(0.0f, 0.0f, 1.0f, 0.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_D)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(-1.0f, 0.0f, 0.0f, 0.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(0.0f, 1.0f, 0.0f, 1.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + (XMVector4Transform(XMVectorSet(0.0f, -1.0f, 0.0f, 0.0f), RotationMatrix) * m_cameraSpeed * static_cast<float>(e.ElapsedTime));
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
    }
    m_CameraYaw -= ImGui::GetMouseDragDelta().x * m_CameraSensibility;
    m_CameraRoll = std::clamp(m_CameraRoll - ImGui::GetMouseDragDelta().y * m_CameraSensibility, (float)-PI / 2.0f, (float)PI / 2.0f);
    ImGui::ResetMouseDragDelta();


    if (!m_freeCamera)
    {
        if (m_autoRotateScene)
        {
            m_autoRotationOffset += static_cast<float>(e.ElapsedTime) * 0.2f;
            if (m_autoRotationOffset > (float)PI * 2.f)
                m_autoRotationOffset -= (float)PI * 2.f;

            m_CameraRoll = lerp(m_CameraRoll, (float)-PI * 0.15f, 2.5f * static_cast<float>(e.ElapsedTime));
            m_CameraYaw = lerp(m_CameraYaw, 0.0f, 2.5f * static_cast<float>(e.ElapsedTime));
        }

        // Compute view matrix
        float4x4 RotationMatrix = XMMatrixRotationY(m_CameraYaw + m_autoRotationOffset) * XMMatrixRotationX(m_CameraRoll);
        float4x4 TranslationMatrix = XMMatrixTranslation(0.f, 0.0f, m_autoCameraDistance);
        m_ViewMatrix = RotationMatrix * TranslationMatrix;
        XMStoreFloat3(&m_cameraPos, XMVector4Transform(XMVectorSet(0.f, 0.0f, m_autoCameraDistance, 1.0f), XMMatrixTranspose(RotationMatrix)));
    }
    else
    {
        // Compute view matrix
        float4x4 RotationMatrix = XMMatrixRotationY(m_CameraYaw + m_autoRotationOffset) * XMMatrixRotationX(m_CameraRoll);
        float4x4 TranslationMatrix = XMMatrixTranslation(m_cameraPos.x, m_cameraPos.y, m_cameraPos.z);
        m_ViewMatrix = TranslationMatrix * RotationMatrix;
    }



    // Update the projection matrix.
    float aspectRatio = GetClientWidth() / static_cast<float>(GetClientHeight());
    m_ProjectionMatrix = XMMatrixPerspectiveFovLH(XMConvertToRadians(m_FoV), aspectRatio, 0.1f, 100000.0f);

    updateImGui();


    
}


// Transition a resource
void ViewDependentMeshletLoD::TransitionResource(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    Microsoft::WRL::ComPtr<ID3D12Resource> resource,
    D3D12_RESOURCE_STATES beforeState, D3D12_RESOURCE_STATES afterState)
{
    CD3DX12_RESOURCE_BARRIER barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        resource.Get(),
        beforeState, afterState);

    commandList->ResourceBarrier(1, &barrier);
}

// Clear a render target.
void ViewDependentMeshletLoD::ClearRTV(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    D3D12_CPU_DESCRIPTOR_HANDLE rtv, FLOAT* clearColor)
{
    commandList->ClearRenderTargetView(rtv, clearColor, 0, nullptr);
}

void ViewDependentMeshletLoD::ClearDepth(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    D3D12_CPU_DESCRIPTOR_HANDLE dsv, FLOAT depth)
{
    commandList->ClearDepthStencilView(dsv, D3D12_CLEAR_FLAG_DEPTH, depth, 0, 0, nullptr);
}

void ViewDependentMeshletLoD::OnRender(RenderEventArgs& e)
{
    super::OnRender(e);

    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_DIRECT);
    auto commandList = commandQueue->GetCommandList();

    UINT currentBackBufferIndex = m_pWindow->GetCurrentBackBufferIndex();
    auto backBuffer = m_pWindow->GetCurrentBackBuffer();
    auto rtv = m_pWindow->GetCurrentRenderTargetView();
    auto dsv = m_DSVHeap->GetCPUDescriptorHandleForHeapStart();


    // Update the MVP matrix
    XMMATRIX mvpMatrix = XMMatrixTranspose(m_ViewMatrix * m_ProjectionMatrix);
    S_Constants constants;

    constants.ViewProjMat = mvpMatrix;
    constants.ViewMat = XMMatrixTranspose(m_ViewMatrix);
    // frustom planes from view-proj-matrix
    // http://gamedevs.org/uploads/fast-extraction-viewing-frustum-planes-from-world-view-projection-matrix.pdf)
    XMFLOAT4X4 m;
    XMStoreFloat4x4(&m, XMMatrixTranspose(mvpMatrix));
    constants.Frustum[0] = float4(m._14 + m._11, m._24 + m._21, m._34 + m._31, m._44 + m._41); // left
    constants.Frustum[1] = float4(m._14 - m._11, m._24 - m._21, m._34 - m._31, m._44 - m._41); // right
    constants.Frustum[2] = float4(m._14 - m._12, m._24 - m._22, m._34 - m._32, m._44 - m._42); // top
    constants.Frustum[3] = float4(m._14 + m._12, m._24 + m._22, m._34 + m._32, m._44 + m._42); // bot
    constants.Frustum[4] = float4(m._13, m._23, m._33, m._43);                                 // near
    constants.Frustum[5] = float4(m._14 - m._13, m._24 - m._23, m._34 - m._33, m._44 - m._43); // far

    constants.CameraWorldPos = m_cameraPos;
    constants.CoTanHalfFoV = m_LoDScale / std::tan(m_FoV / 2.0f);
    constants.CurrTime = (float)m_totalRunTime;
    constants.shadingSelection = m_shadingMode;
    constants.BoolConstants = 0;
    constants.SceneObjectCount = (uint)m_scene.m_scene_objects.size();
    if (m_frustumCulling) constants.BoolConstants |= FRUSTUM_CULLING_BIT_POS;
    
    memcpy(m_mappedConstantData, &constants, sizeof(S_Constants));

    // Clear the render targets.
    {
        TransitionResource(commandList, backBuffer,
            D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET);

        ClearRTV(commandList, rtv, m_ClearColor);
        ClearDepth(commandList, dsv);
    }

    // main render pass
    commandList->SetPipelineState(m_PipelineState.Get());
    commandList->SetGraphicsRootSignature(m_RootSignature.Get());

    commandList->RSSetViewports(1, &m_Viewport);
    commandList->RSSetScissorRects(1, &m_ScissorRect);

    commandList->OMSetRenderTargets(1, &rtv, FALSE, &dsv);

    ID3D12DescriptorHeap* heaps[] = { m_CBV_SRV_UAV_Heap.Get() };
    commandList->SetDescriptorHeaps(_countof(heaps), heaps);

    //set constants
    commandList->SetGraphicsRootConstantBufferView(0, m_ConstantsBuffer->GetGPUVirtualAddress());
    // set bindless vertex buffers
    commandList->SetGraphicsRootDescriptorTable(1, m_VerticesSrvHandle);
    // set bindless vertex-indices buffers
    commandList->SetGraphicsRootDescriptorTable(2, m_VertexIndicesSrvHandle);
    // set bindless primitive-indices buffers
    commandList->SetGraphicsRootDescriptorTable(3, m_PrimitiveIndicesSrvHandle);
    // set bindless meshlet buffers
    commandList->SetGraphicsRootDescriptorTable(4, m_MeshletsSrvHandle);
    // set bindless meshlet groups buffers
    commandList->SetGraphicsRootDescriptorTable(5, m_MeshletGroupsSrvHandle);
    // set scene objects buffer
    commandList->SetGraphicsRootShaderResourceView(6, m_ObjectsBuffer.Get()->GetGPUVirtualAddress());
    // set work queue buffer
    commandList->SetGraphicsRootUnorderedAccessView(7, m_WorkQueueBuffer.Get()->GetGPUVirtualAddress());
    // clear with zeros and set work queue counters buffer
    commandList->CopyResource(m_WorkQueueCountersBuffer.Get(), m_WorkQueueCountersClearValuesBuffer.Get());
    commandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::UAV(m_WorkQueueCountersBuffer.Get()));
    commandList->SetGraphicsRootUnorderedAccessView(8, m_WorkQueueCountersBuffer.Get()->GetGPUVirtualAddress());
    
    
    
    commandList->DispatchMesh(PERSISTENT_THREAD_COUNT, 1, 1);

    // Render Dear ImGui graphics
    commandList->OMSetRenderTargets(1, &rtv, FALSE, nullptr);
    commandList->SetDescriptorHeaps(1, m_ImGuiDescriptorHeap.GetAddressOf());
    ImGui_ImplDX12_RenderDrawData(ImGui::GetDrawData(), commandList.Get());
 

    // Present
    {
        TransitionResource(commandList, backBuffer,
            D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PRESENT);

        m_FenceValues[currentBackBufferIndex] = commandQueue->ExecuteCommandList(commandList);
        
        currentBackBufferIndex = m_pWindow->Present();

        commandQueue->WaitForFenceValue(m_FenceValues[currentBackBufferIndex]);
    }
}

void ViewDependentMeshletLoD::OnKeyPressed(KeyEventArgs& e)
{
    super::OnKeyPressed(e);

    switch (e.Key)
    {
    case KeyCode::Escape:
        Application::Get().Quit(0);
        break;
    case KeyCode::Enter:
        if (e.Alt)
        {
    case KeyCode::F11:
        m_pWindow->ToggleFullscreen();
        break;
        }
    case KeyCode::V:
        m_pWindow->ToggleVSync();
        break;
    }
}

void ViewDependentMeshletLoD::OnMouseWheel(MouseWheelEventArgs& e)
{
    m_autoCameraDistance = std::clamp(m_autoCameraDistance - e.WheelDelta * m_CameraScrollScale, 1.0f, 200.0f);
}

void ViewDependentMeshletLoD::initImGui()
{
    static bool alreadyInitialized = false;
    if (alreadyInitialized) return;

    // Simple free list based allocator
    struct ExampleDescriptorHeapAllocator
    {
        ID3D12DescriptorHeap* Heap = nullptr;
        D3D12_DESCRIPTOR_HEAP_TYPE  HeapType = D3D12_DESCRIPTOR_HEAP_TYPE_NUM_TYPES;
        D3D12_CPU_DESCRIPTOR_HANDLE HeapStartCpu;
        D3D12_GPU_DESCRIPTOR_HANDLE HeapStartGpu;
        UINT                        HeapHandleIncrement;
        ImVector<int>               FreeIndices;

        void Create(ID3D12Device* device, ID3D12DescriptorHeap* heap)
        {
            IM_ASSERT(Heap == nullptr && FreeIndices.empty());
            Heap = heap;
            D3D12_DESCRIPTOR_HEAP_DESC desc = heap->GetDesc();
            HeapType = desc.Type;
            HeapStartCpu = Heap->GetCPUDescriptorHandleForHeapStart();
            HeapStartGpu = Heap->GetGPUDescriptorHandleForHeapStart();
            HeapHandleIncrement = device->GetDescriptorHandleIncrementSize(HeapType);
            FreeIndices.reserve((int)desc.NumDescriptors);
            for (int n = desc.NumDescriptors; n > 0; n--)
                FreeIndices.push_back(n);
        }
        void Destroy()
        {
            Heap = nullptr;
            FreeIndices.clear();
        }
        void Alloc(D3D12_CPU_DESCRIPTOR_HANDLE* out_cpu_desc_handle, D3D12_GPU_DESCRIPTOR_HANDLE* out_gpu_desc_handle)
        {
            IM_ASSERT(FreeIndices.Size > 0);
            int idx = FreeIndices.back();
            FreeIndices.pop_back();
            out_cpu_desc_handle->ptr = HeapStartCpu.ptr + (idx * HeapHandleIncrement);
            out_gpu_desc_handle->ptr = HeapStartGpu.ptr + (idx * HeapHandleIncrement);
        }
        void Free(D3D12_CPU_DESCRIPTOR_HANDLE out_cpu_desc_handle, D3D12_GPU_DESCRIPTOR_HANDLE out_gpu_desc_handle)
        {
            int cpu_idx = (int)((out_cpu_desc_handle.ptr - HeapStartCpu.ptr) / HeapHandleIncrement);
            int gpu_idx = (int)((out_gpu_desc_handle.ptr - HeapStartGpu.ptr) / HeapHandleIncrement);
            IM_ASSERT(cpu_idx == gpu_idx);
            FreeIndices.push_back(cpu_idx);
        }
    };





    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();


    

    D3D12_DESCRIPTOR_HEAP_DESC descHeap = {};
    descHeap.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    descHeap.NumDescriptors = 64;
    descHeap.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    descHeap.NodeMask = 0;

    ThrowIfFailed(Application::Get().GetDevice().Get()->CreateDescriptorHeap(
        &descHeap, IID_PPV_ARGS(&m_ImGuiDescriptorHeap)
    ));

    
    // Setup Platform/Renderer backends
    ImGui_ImplWin32_Init(m_pWindow->GetWindowHandle());

    static std::shared_ptr<CommandQueue> commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_DIRECT);
    static ExampleDescriptorHeapAllocator g_pd3dSrvDescHeapAlloc;
    g_pd3dSrvDescHeapAlloc.Create(Application::Get().GetDevice().Get(), m_ImGuiDescriptorHeap.Get());
    

    ImGui_ImplDX12_InitInfo init_info = {};
    init_info.Device = Application::Get().GetDevice().Get();
    init_info.CommandQueue = commandQueue->GetD3D12CommandQueue().Get();
    init_info.NumFramesInFlight = 3;
    init_info.RTVFormat = DXGI_FORMAT_R8G8B8A8_UNORM;
    init_info.DSVFormat = DXGI_FORMAT_UNKNOWN;
    // Allocating SRV descriptors (for textures) is up to the application, so we provide callbacks.
    // (current version of the backend will only allocate one descriptor, future versions will need to allocate more)
    init_info.SrvDescriptorHeap = m_ImGuiDescriptorHeap.Get();
    init_info.SrvDescriptorAllocFn = [](ImGui_ImplDX12_InitInfo*, D3D12_CPU_DESCRIPTOR_HANDLE* out_cpu_handle, D3D12_GPU_DESCRIPTOR_HANDLE* out_gpu_handle) { return g_pd3dSrvDescHeapAlloc.Alloc(out_cpu_handle, out_gpu_handle); };
    init_info.SrvDescriptorFreeFn = [](ImGui_ImplDX12_InitInfo*, D3D12_CPU_DESCRIPTOR_HANDLE cpu_handle, D3D12_GPU_DESCRIPTOR_HANDLE gpu_handle) { return g_pd3dSrvDescHeapAlloc.Free(cpu_handle, gpu_handle); };
    ImGui_ImplDX12_Init(&init_info);


    alreadyInitialized = true;
}


void ViewDependentMeshletLoD::updateImGui()
{
    
    DXGI_QUERY_VIDEO_MEMORY_INFO videoMemoryInfo;
    Application::Get().m_dxgiAdapter->QueryVideoMemoryInfo(0, DXGI_MEMORY_SEGMENT_GROUP_LOCAL, &videoMemoryInfo);


    PROCESS_MEMORY_COUNTERS_EX memoryInfo;
    memoryInfo.cb = sizeof(memoryInfo);

    //MEMORYSTATUSEX memoryInfo;
    //memoryInfo.dwLength = sizeof(MEMORYSTATUSEX);

    const unsigned int fpsHistorySize = 250;
    static float fpsHistory[fpsHistorySize] = { 0 };
    static float fpsHistoryOffset = 0.0f;
    float fpsHistorySpeed = 50.0f;

    fpsHistory[(unsigned int)fpsHistoryOffset] = std::max(fpsHistory[(unsigned int)fpsHistoryOffset], (float)m_frameTime);
    fpsHistory[((unsigned int)fpsHistoryOffset + 1) % fpsHistorySize] = 0;
    fpsHistoryOffset += clamp(fpsHistorySpeed * (float)m_frameTime, 0.0f, 1.0f);
    if (fpsHistoryOffset >= fpsHistorySize) fpsHistoryOffset = 0.0f;



    // Start the Dear ImGui frame
    
    ImGui_ImplWin32_NewFrame();
    ImGui_ImplDX12_NewFrame();
    ImGui::NewFrame();


    ImGui::SetNextWindowSize(ImVec2(375, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 375 - 12), 12.0f), ImGuiCond_Always);
    ImGui::Begin("Settings");                           
   

    if (!ImGui::CollapsingHeader("Render Settings"))
    {
        ImGui::Checkbox("Meshlet based Frustum Culling", &m_frustumCulling);
        ImGui::Checkbox("Object Culling", &m_objectCulling);
        ImGui::Checkbox("Level Of Detail", &m_LoD_Enabled);
        ImGui::SliderFloat("LoD Scale", &m_LoDScale, 0.01f, 10.0f);
        ImGui::ColorEdit4("Clear Color", m_ClearColor);
        if (ImGui::Button("Toggle Fullscreen"))
        {
            m_pWindow->ToggleFullscreen();
        }
        if (ImGui::Button("Toggle V-Sync"))
        {
            m_pWindow->ToggleVSync();
        }
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Camera Settings"))
    {
        ImGui::Checkbox("Free Camera", &m_freeCamera);
        if (m_freeCamera)
        {
            ImGui::SliderFloat("Camera Speed", &m_cameraSpeed, 0.1f, 100.0f);
        }
        else
        {
            ImGui::SliderFloat("Scroll Scale", &m_CameraScrollScale, 0.1f, 5.0f);
            ImGui::Checkbox("Auto Rotation", &m_autoRotateScene);
        }
        ImGui::SliderFloat("FoV", &m_FoV, 45.0f, 120.0f);
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Model Settings"))
    {
        ImGui::InputText("", m_model_file_path, sizeof(m_model_file_path));
        ImGui::SameLine();
        if (ImGui::Button("Load Model"))
        {
            Application::Get().Flush();
            UnloadContent();
            m_scene.init(m_model_file_path);
            LoadContent();
        }
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Debug Settings"))
    {
        if (ImGui::Checkbox("Wireframe", &m_wireframe)) createPSO();
        ImGui::SameLine();
        if (ImGui::Checkbox("Back Face Culling", &m_backFaceCulling)) createPSO();
        
        
        static int selected = 0;  // Index of the selected option
        const char* options[] = { "Disable Debug Visals", "Show Meshlets", "Show LoDs"};
        for (int i = 0; i < IM_ARRAYSIZE(options); i++) {
            if (ImGui::RadioButton(options[i], selected == i)) {
                selected = i;  // Update selection when clicked
                m_shadingMode = static_cast<ShadingMode>(i);
            }
        }
        
        if (ImGui::Button("Force Lag"))
            Sleep(1000);
        
        
    }
    ImGui::Spacing();

    ImGui::End();

    // Performance Statistics
    ImGui::SetNextWindowSize(ImVec2(260, 0));
    ImGui::SetNextWindowPos(ImVec2(12.0f, 12.0f), ImGuiCond_Always);
    ImGui::Begin("Performance Statistics");

    if (!ImGui::CollapsingHeader("Run-Time Perfromance Statistics"))
    {
        ImGui::Text("%.1f FPS", m_fps);
        ImGui::Text("%.2fms  ", m_frameTime * 1000);
        ImGui::SameLine();
        ImGui::PlotLines("##IDfix", fpsHistory, fpsHistorySize);
        if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memoryInfo, sizeof(memoryInfo)))
            ImGui::Text("%.3fGB RAM Usage", (memoryInfo.WorkingSetSize) / 1073741824.0f);
        ImGui::Text("%.3fGB VRAM Usage", videoMemoryInfo.CurrentUsage / 1000000000.0f);
    }
    if (!ImGui::CollapsingHeader("Scene Stats"))
    {
        ImGui::Text("scene  objects count: %d", (unsigned int)m_scene.m_scene_objects.size());
        ImGui::Text("unique  meshes count: %d", (unsigned int)m_scene.m_mesh_count);


        if (!ImGui::CollapsingHeader("Scene Processing Times"))
        {
            ImGui::Text("model load time:         %.4f sec", 0);
            ImGui::Text("LoD generation time:     %.4f sec", 0);
            ImGui::Text("meshlet generation time: %.4f sec", 0);
            ImGui::Separator();
            ImGui::Text("total load time:         %.4f sec", 0);
        }

    }

    ImGui::End();


    // Help Window
    ImGui::SetNextWindowSize(ImVec2(275, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 275 - 375 - 12 - 7), 12.0f), ImGuiCond_Always);
    ImGui::Begin("Help");

    if (!ImGui::CollapsingHeader("Controls"))
    {
        ImGui::Text("W,A,S,D - move camera horizontally");
        ImGui::Text("shift   - move camera up");
        ImGui::Text("ctrl    - move camera down");

        ImGui::Text("esc     - close application");
        ImGui::Text("F11     - toggle fullscreen");
        ImGui::Text("V       - toggle v-sync");
    } 
    ImGui::End();

    ImGui::Render(); 
}