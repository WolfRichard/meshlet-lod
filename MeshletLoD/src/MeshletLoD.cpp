#include <directx/d3dx12.h>
#include <MeshletLoD.h>

#include <psapi.h>

#include "HLSLnames.h"

#include "../shaders/Structures.fxh"

#include <Application.h>
#include <CommandQueue.h>
#include <Helpers.h>
#include <Window.h>
#if defined(min)
#undef min
#endif

#if defined(max)
#undef max
#endif
#include <sstream>

#include <imgui.h>
#include "imgui_impl_win32.h"
#include "imgui_impl_dx12.h"
#include "imgui_internal.h"

#include <wrl.h>
using namespace Microsoft::WRL;

#include <d3dcompiler.h>
#include <algorithm> 
#include <DirectXTex.h>



using namespace DirectX;


// Clamp a value between a min and max range.
template<typename T>
constexpr const T& clamp(const T& val, const T& min, const T& max)
{
    return val < min ? min : val > max ? max : val;
}

MeshletLoD::MeshletLoD(const std::wstring& name, int width, int height, bool vSync)
    : super(name, width, height, vSync)
    , m_ScissorRect(CD3DX12_RECT(0, 0, LONG_MAX, LONG_MAX))
    , m_Viewport(CD3DX12_VIEWPORT(0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height)))
    , m_FoV(90.0)
    , m_ContentLoaded(false)
{
    m_scene.init(m_model_file_path);
}

void MeshletLoD::UpdateBufferResource(
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

void MeshletLoD::createPSO()
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


void MeshletLoD::setupConstantsUploadBuffer()
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
void MeshletLoD::setupSrvAndBuffer(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
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
void MeshletLoD::setupBindlessSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle, 
                                                         D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle, 
                                                         D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                                         std::vector<std::vector<T>>& cpuBuffers, 
                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                                         std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                                         Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                                         unsigned int descriptorSize)
{
    auto device = Application::Get().GetDevice();
    srvGpuHandle = nextAvailableGpuSrvHandle;
    
    for (uint i = 0; i < cpuBuffers.size(); i++)
    {
        std::vector<T>& cpuSingleBuffer = cpuBuffers[i];

        gpuBuffers.push_back(ComPtr<ID3D12Resource>());
        copyBuffers.push_back(ComPtr<ID3D12Resource>());
        
        UpdateBufferResource(commandList, gpuBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)cpuSingleBuffer.size(), sizeof(T), cpuSingleBuffer.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)cpuSingleBuffer.size();
        srvDesc.Buffer.StructureByteStride = sizeof(T);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(gpuBuffers.back().Get(), &srvDesc, nextAvailableCpuSrvHandle);
        nextAvailableCpuSrvHandle.ptr += descriptorSize;
        nextAvailableGpuSrvHandle.ptr += descriptorSize;
    }
}

// compatible with tight packing of 4 chars into 1 unsigned integers
void MeshletLoD::setupBindlessUCharToUIntSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                                                                    D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle,
                                                                    D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                                                    std::vector<std::vector<unsigned char>>& unsignedCharsCpuBuffers,
                                                                    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                                                    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                                                    Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                                                    unsigned int descriptorSize)
{
    auto device = Application::Get().GetDevice();
    srvGpuHandle = nextAvailableGpuSrvHandle;

    for (uint i = 0; i < unsignedCharsCpuBuffers.size(); i++)
    {
        std::vector<unsigned char>& cpuSingleBuffer = unsignedCharsCpuBuffers[i];
        // fill with zeros to allign its size to 32bits
        while (cpuSingleBuffer.size() % 4 != 0)
            cpuSingleBuffer.push_back(0);

        gpuBuffers.push_back(ComPtr<ID3D12Resource>());
        copyBuffers.push_back(ComPtr<ID3D12Resource>());

        UpdateBufferResource(commandList, gpuBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)cpuSingleBuffer.size() / 4, sizeof(uint), cpuSingleBuffer.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)cpuSingleBuffer.size() / 4;
        srvDesc.Buffer.StructureByteStride = sizeof(uint);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(gpuBuffers.back().Get(), &srvDesc, nextAvailableCpuSrvHandle);
        nextAvailableCpuSrvHandle.ptr += descriptorSize;
        nextAvailableGpuSrvHandle.ptr += descriptorSize;
    }
}



bool MeshletLoD::LoadContent()
{
    auto api_setup_time_start = std::chrono::high_resolution_clock::now();

    auto device = Application::Get().GetDevice();
    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_COPY);
    auto commandList = commandQueue->GetCommandList();

    setupConstantsUploadBuffer();

    // height map texture sampler
    D3D12_DESCRIPTOR_HEAP_DESC samplerHeapDesc = {};
    samplerHeapDesc.NumDescriptors = 1;
    samplerHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_SAMPLER;
    samplerHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    device->CreateDescriptorHeap(&samplerHeapDesc, IID_PPV_ARGS(&m_Sampler_Heap));

    D3D12_SAMPLER_DESC samplerDesc = {};
    samplerDesc.Filter = D3D12_FILTER_MIN_MAG_MIP_LINEAR;
    samplerDesc.AddressU = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
    samplerDesc.AddressV = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
    samplerDesc.AddressW = D3D12_TEXTURE_ADDRESS_MODE_WRAP;
    samplerDesc.MinLOD = 0;
    samplerDesc.MaxLOD = D3D12_FLOAT32_MAX;
    samplerDesc.MaxAnisotropy = 1;

    // Create the sampler
    device->CreateSampler(&samplerDesc, m_Sampler_Heap->GetCPUDescriptorHandleForHeapStart());

    // setup cbv/srv/uav descriptor heap
    D3D12_DESCRIPTOR_HEAP_DESC gpuHeapDesc = {};
    gpuHeapDesc.NumDescriptors = 1                           // single scene objects buffer
                               + 1                           // single work queue buffer
                               + 1                           // single work queue counter buffer
                               + 1                           // single work queue counter clear values buffer
                               + 1                           // global mesh payload buffer
                               + 1                           // height map texture
                               + (uint)m_scene.m_mesh_count  // meshlets buffer per unique mesh
                               + (uint)m_scene.m_mesh_count  // morph indices buffer per unique mesh
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


    // load height map texture
    ScratchImage height_map_test_image;
    HRESULT hr = LoadFromWICFile(L"./assets/textures/heightmap.png", WIC_FLAGS_NONE, nullptr, height_map_test_image);
    if (FAILED(hr)) {
        OutputDebugString("FAILED TO LOAD TEST TEXTURE!\n");
        assert(false);
    }
    const Image* img = height_map_test_image.GetImage(0, 0, 0);

    D3D12_RESOURCE_DESC textureDesc = {};
    textureDesc.MipLevels = 0;
    textureDesc.Format = img->format;
    textureDesc.Width = (uint)img->width;
    textureDesc.Height = (uint)img->height;
    textureDesc.Flags = D3D12_RESOURCE_FLAG_NONE;
    textureDesc.DepthOrArraySize = 1;
    textureDesc.SampleDesc.Count = 1;
    textureDesc.SampleDesc.Quality = 0;
    textureDesc.Dimension = D3D12_RESOURCE_DIMENSION_TEXTURE2D;

    copyBuffers.push_back(ComPtr<ID3D12Resource>());

    UINT64 bufferSize = 0;
    UINT numRows = 0;
    UINT64 rowSizeInBytes = 0;
    D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint = {};
    device->GetCopyableFootprints(&textureDesc, 0, 1, 0, &footprint,
        &numRows, &rowSizeInBytes, &bufferSize
    );
 
    // Create a committed resource for the GPU resource in a default heap.
    ThrowIfFailed(device->CreateCommittedResource(
        &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_DEFAULT),
        D3D12_HEAP_FLAG_NONE,
        &textureDesc,
        D3D12_RESOURCE_STATE_COMMON,
        nullptr,
        IID_PPV_ARGS(m_HeightMapTexture.GetAddressOf())));

    // Create an committed resource for the upload.
    if (img->pixels)
    {
        ThrowIfFailed(device->CreateCommittedResource(
            &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
            D3D12_HEAP_FLAG_NONE,
            &CD3DX12_RESOURCE_DESC::Buffer(bufferSize),
            D3D12_RESOURCE_STATE_GENERIC_READ,
            nullptr,
            IID_PPV_ARGS(copyBuffers.back().GetAddressOf())));

        D3D12_SUBRESOURCE_DATA subresourceData = {};
        subresourceData.pData = img->pixels;
        subresourceData.RowPitch = img->rowPitch;
        subresourceData.SlicePitch = img->slicePitch;

        UpdateSubresources(commandList.Get(),
            *m_HeightMapTexture.GetAddressOf(), *copyBuffers.back().GetAddressOf(),
            0, 0, 1, &subresourceData);
    }
    else
        assert(false); // currently no backup texture if height-map loading fails

    D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
    srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    srvDesc.Format = textureDesc.Format;
    srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
    srvDesc.Texture2D.MostDetailedMip = 0;
    srvDesc.Texture2D.MipLevels = height_map_test_image.GetMetadata().mipLevels;

    device->CreateShaderResourceView(m_HeightMapTexture.Get(), &srvDesc, nextCpuSrvHandle);
    m_HeightMapTextureSrvHandle = nextGpuSrvHandle;

    nextCpuSrvHandle.ptr += descriptorSize;
    nextGpuSrvHandle.ptr += descriptorSize;
    
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
    setupBindlessUCharToUIntSrvAndBuffers(m_PrimitiveIndicesSrvHandle,
                                          nextGpuSrvHandle, nextCpuSrvHandle,
                                          m_scene.m_primitive_indices, m_PrimitiveIndicesBuffers,
                                          copyBuffers, commandList, descriptorSize);

    // morph-indices buffers
    setupBindlessSrvAndBuffers(m_MorphIndicesSrvHandle,
                               nextGpuSrvHandle, nextCpuSrvHandle,
                               m_scene.m_morph_indices, m_MorphIndicesBuffers,
                               copyBuffers, commandList, descriptorSize);

    // meshlets buffers
    setupBindlessSrvAndBuffers(m_MeshletsSrvHandle,
                               nextGpuSrvHandle, nextCpuSrvHandle,
                               m_scene.m_meshlets, m_MeshletBuffers,
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
    std::vector<S_WorkQueueCounters> temporaryCpuWorkQueueCountersData(1, { 0, 0, 0 });
    setupSrvAndBuffer(m_WorkQueueCountersSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      temporaryCpuWorkQueueCountersData, m_WorkQueueCountersBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);
    setupSrvAndBuffer(m_WorkQueueCountersClearValuesSrvHandle,
                      nextGpuSrvHandle, nextCpuSrvHandle,
                      temporaryCpuWorkQueueCountersData, m_WorkQueueCountersClearValuesBuffer,
                      copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_NONE);

    // work queue counters buffer and clear values
    std::vector<S_PayloadEntry> temporaryMeshPaloadData(MAX_DISPATCH_MESH_GROUP_COUNT, { 0, 0, 0, 0 });
    setupSrvAndBuffer(m_GlobalMeshPayloadSrvHandle,
        nextGpuSrvHandle, nextCpuSrvHandle,
        temporaryMeshPaloadData, m_GlobalMeshPayloadBuffer,
        copyBuffers, commandList, descriptorSize, D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);

    // Create the descriptor heap for the depth-stencil view.
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&m_DSVHeap)));

    // Load shaders.
    ThrowIfFailed(D3DReadFileToBlob(L"PixelShader.cso", &m_pixelShaderBlob));
    ThrowIfFailed(D3DReadFileToBlob(L"MeshShader.cso", &m_meshShaderBlob));
    ThrowIfFailed(D3DReadFileToBlob(L"TaskShader.cso", &m_taskShaderBlob));

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

    CD3DX12_ROOT_PARAMETER1 rootParameters[12];

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

    // morph indices
    CD3DX12_DESCRIPTOR_RANGE1 srvMorphIndicesRange;
    srvMorphIndicesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_mesh_count, 0, 5);
    rootParameters[5].InitAsDescriptorTable(1, &srvMorphIndicesRange, D3D12_SHADER_VISIBILITY_MESH);

    // scene objects
    rootParameters[6].InitAsShaderResourceView(0, 0);

    // work queue objects
    rootParameters[7].InitAsUnorderedAccessView(0, 0);

    // work queue counters objects
    rootParameters[8].InitAsUnorderedAccessView(1, 0);

    // global mesh payload
    rootParameters[9].InitAsUnorderedAccessView(2, 0);
     
    // height map texture
    CD3DX12_DESCRIPTOR_RANGE1 textureRange;
    textureRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 1, 0);

    rootParameters[10].InitAsDescriptorTable(1, &textureRange, D3D12_SHADER_VISIBILITY_ALL);

    // height map texture sampler
    CD3DX12_DESCRIPTOR_RANGE1 samplerRange;
    samplerRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SAMPLER, 1, 0, 0);
    rootParameters[11].InitAsDescriptorTable(1, &samplerRange, D3D12_SHADER_VISIBILITY_ALL);

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


    auto api_setup_time_end = std::chrono::high_resolution_clock::now();
    m_apiSetupTime = api_setup_time_end - api_setup_time_start;

    return true;
}

void MeshletLoD::ResizeDepthBuffer(int width, int height)
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

void MeshletLoD::OnResize(ResizeEventArgs& e)
{
    if (e.Width != GetClientWidth() || e.Height != GetClientHeight())
    {
        super::OnResize(e);

        m_Viewport = CD3DX12_VIEWPORT(0.0f, 0.0f,
            static_cast<float>(e.Width), static_cast<float>(e.Height));

        ResizeDepthBuffer(e.Width, e.Height);
    }
}

void MeshletLoD::UnloadContent()
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
    m_MorphIndicesBuffers.clear();
    m_ObjectsBuffer.Reset();
    m_ConstantsBuffer.Reset();
    m_WorkQueueBuffer.Reset();
    m_WorkQueueCountersBuffer.Reset();
    m_GlobalMeshPayloadBuffer.Reset();

    m_commandSignature.Reset();
}

void MeshletLoD::OnUpdate(UpdateEventArgs& e)
{
    // update frame time and fps
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
        char buffer[512];
        sprintf_s(buffer, "FPS: %f\n", m_fps);
        OutputDebugStringA(buffer);

        frameCount = 0;
        totalTime = 0.0;
    }
    
    // Camera Controll
    if (m_freeCamera)
    {
        // Compute relative directions from yaw and pitch
        XMVECTOR forward = XMVectorSet(cosf(m_CameraRoll) * sinf(m_CameraYaw - m_autoRotationOffset), sinf(m_CameraRoll), cosf(m_CameraRoll) * cosf(m_CameraYaw - m_autoRotationOffset), 0.0f);
        XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);  
        XMVECTOR right = XMVector3Normalize(XMVector3Cross(up, forward));
        up = XMVector3Normalize(XMVector3Cross(forward, right));

        if (ImGui::IsKeyDown(ImGuiKey_W)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + forward * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_A)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) - right * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_S)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) - forward * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_D)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + right * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) - up * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
        if (ImGui::IsKeyDown(ImGuiKey_LeftShift)) {
            XMVECTOR newCamPos = XMLoadFloat3(&m_cameraPos) + up * m_cameraSpeed * static_cast<float>(e.ElapsedTime);
            XMStoreFloat3(&m_cameraPos, newCamPos);
        }
    }
    m_CameraYaw += ImGui::GetMouseDragDelta().x * m_CameraSensibility;
    m_CameraRoll = std::clamp(m_CameraRoll - ImGui::GetMouseDragDelta().y * m_CameraSensibility, (float)-PI / 2.001f, (float)PI / 2.001f);
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
        XMVECTOR forward = XMVectorSet(cosf(m_CameraRoll) * sinf(m_CameraYaw - m_autoRotationOffset), sinf(m_CameraRoll), cosf(m_CameraRoll) * cosf(m_CameraYaw - m_autoRotationOffset), 0.0f);
        XMStoreFloat3(&m_cameraPos, -forward * m_autoCameraDistance);
        m_ViewMatrix = XMMatrixLookAtLH(XMLoadFloat3(&m_cameraPos), XMVectorSet(0, 0, 0, 0), XMVectorSet(0, 1, 0, 0));
    }
    else
    {
        // Compute forward direction from yaw and pitch
        XMVECTOR forward = XMVectorSet(cosf(m_CameraRoll) * sinf(m_CameraYaw - m_autoRotationOffset), sinf(m_CameraRoll), cosf(m_CameraRoll) * cosf(m_CameraYaw - m_autoRotationOffset), 0.0f);
        // Normalize the forward vector
        forward = XMVector3Normalize(forward);
        // Load camera position into XMVECTOR
        XMVECTOR eye = XMLoadFloat3(&m_cameraPos);
        // Compute focus point (camera position + forward * distance)
        XMVECTOR focus = XMVectorAdd(eye, forward);
        m_ViewMatrix = XMMatrixLookAtLH(XMLoadFloat3(&m_cameraPos), focus, XMVectorSet(0, 1, 0, 0));
    }

    // Update the projection matrix.
    float aspectRatio = GetClientWidth() / static_cast<float>(GetClientHeight());
    m_ProjectionMatrix = XMMatrixPerspectiveFovLH(XMConvertToRadians(m_FoV), aspectRatio, 0.001f, 100000.0f);

    updateImGui();
}


// Transition a resource
void MeshletLoD::TransitionResource(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    Microsoft::WRL::ComPtr<ID3D12Resource> resource,
    D3D12_RESOURCE_STATES beforeState, D3D12_RESOURCE_STATES afterState)
{
    CD3DX12_RESOURCE_BARRIER barrier = CD3DX12_RESOURCE_BARRIER::Transition(
        resource.Get(),
        beforeState, afterState);

    commandList->ResourceBarrier(1, &barrier);
}

void MeshletLoD::ClearRTV(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    D3D12_CPU_DESCRIPTOR_HANDLE rtv, FLOAT* clearColor)
{
    commandList->ClearRenderTargetView(rtv, clearColor, 0, nullptr);
}

void MeshletLoD::ClearDepth(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
    D3D12_CPU_DESCRIPTOR_HANDLE dsv, FLOAT depth)
{
    commandList->ClearDepthStencilView(dsv, D3D12_CLEAR_FLAG_DEPTH, depth, 0, 0, nullptr);
}

void MeshletLoD::OnRender(RenderEventArgs& e)
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
    constants.ProjMat = XMMatrixTranspose(m_ProjectionMatrix);
    constants.ScreenWidth = GetClientWidth();
    constants.ScreenHeight = GetClientHeight();

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

    float scale_x = sqrtf(m._11 * m._11 + m._21 * m._21 + m._31 * m._31);
    float scale_y = sqrtf(m._12 * m._12 + m._22 * m._22 + m._32 * m._32);
    float scale_z = sqrtf(m._13 * m._13 + m._23 * m._23 + m._33 * m._33);
    constants.MaxScaleFactor_ViewProjMat = std::max(scale_x, std::max(scale_y, scale_z));

    // do not update camera position constant when camera is locked for debug purposes
    if (m_lockCameraShaderConstant)
        constants.CameraWorldPos = m_lockedCameraPos;
    else
        constants.CameraWorldPos = m_cameraPos;
    
    constants.CoTanHalfFoV = 1 / std::tan(m_FoV / 2.0f);
    constants.LoD_Scale = m_LoDScale;
    constants.CurrTime = (float)m_totalRunTime;
    constants.shadingSelection = m_shadingMode;
    constants.BoolConstants = 0;
    constants.SceneObjectCount = (uint)m_scene.m_scene_objects.size();
    constants.DebugFloatSliderValue = m_debugFloatSlider;
    constants.TriPlanarMappingScale = m_triplanarScale;
    constants.TriPlanarBlendGrade = m_triplanarBlendGrade;
    constants.HeightMapDisplacementScale = m_displacementScale;

    //bool constants
    if (m_frustumCulling)           constants.BoolConstants |= FRUSTUM_CULLING_BIT_POS;
    if (m_geo_morphing)             constants.BoolConstants |= GEO_MORPHING_BIT_POS;
    if (m_screen_space_LoD)         constants.BoolConstants |= SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS;
    if (m_tessellation)             constants.BoolConstants |= TRESSELLATION_BIT_POS;
    if (m_triplanarMapping)         constants.BoolConstants |= TRI_PLANAR_TEXTURE_MAPPING_BIT_POS;
    if (m_allowLighting)            constants.BoolConstants |= NORMAL_LIGHTING_BIT_POS;

    // push new constants data from cpu to gpu
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

    ID3D12DescriptorHeap* heaps[] = { m_CBV_SRV_UAV_Heap.Get(), m_Sampler_Heap.Get()};
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
    // set bindless morph indices buffers
    commandList->SetGraphicsRootDescriptorTable(5, m_MorphIndicesSrvHandle);
    // set scene objects buffer
    commandList->SetGraphicsRootShaderResourceView(6, m_ObjectsBuffer.Get()->GetGPUVirtualAddress());
    // set work queue buffer
    commandList->SetGraphicsRootUnorderedAccessView(7, m_WorkQueueBuffer.Get()->GetGPUVirtualAddress());
    // clear with zeros and set work queue counters buffer
    commandList->CopyResource(m_WorkQueueCountersBuffer.Get(), m_WorkQueueCountersClearValuesBuffer.Get());
    commandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::UAV(m_WorkQueueCountersBuffer.Get()));
    commandList->SetGraphicsRootUnorderedAccessView(8, m_WorkQueueCountersBuffer.Get()->GetGPUVirtualAddress());
    // set global mesh payload buffer
    commandList->SetGraphicsRootUnorderedAccessView(9, m_GlobalMeshPayloadBuffer.Get()->GetGPUVirtualAddress());
    // set height map texture
    commandList->SetGraphicsRootDescriptorTable(10,m_HeightMapTextureSrvHandle);
    //commandList->SetGraphicsRootShaderResourceView(10, m_HeightMapTexture.Get()->GetGPUVirtualAddress());
    // set hight map texture sampler
    commandList->SetGraphicsRootDescriptorTable(11, m_Sampler_Heap->GetGPUDescriptorHandleForHeapStart());

    
    // dispatch a fixed amount of persistent threads to process and render the scene
    assert((PERSISTENT_THREAD_COUNT % GROUP_SIZE) == 0);
    commandList->DispatchMesh(PERSISTENT_THREAD_COUNT / GROUP_SIZE, 1, 1);


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

void MeshletLoD::OnKeyPressed(KeyEventArgs& e)
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

void MeshletLoD::OnMouseWheel(MouseWheelEventArgs& e)
{
    m_autoCameraDistance = std::clamp(m_autoCameraDistance - e.WheelDelta * m_CameraScrollScale, 1.0f, 200.0f);
}

void MeshletLoD::initImGui()
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
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

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


void MeshletLoD::updateImGui()
{
    // update VRAM usage
    DXGI_QUERY_VIDEO_MEMORY_INFO videoMemoryInfo;
    Application::Get().m_dxgiAdapter->QueryVideoMemoryInfo(0, DXGI_MEMORY_SEGMENT_GROUP_LOCAL, &videoMemoryInfo);

    // update RAM usage
    PROCESS_MEMORY_COUNTERS_EX memoryInfo;
    memoryInfo.cb = sizeof(memoryInfo);

    // keep track of recent frame time
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

    // render and debug settings window
    ImGui::SetNextWindowSize(ImVec2(375, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 375 - 12), 12.0f), ImGuiCond_Always);
    ImGui::Begin("Settings");                           
    if (!ImGui::CollapsingHeader("Render Settings"))
    {
        ImGui::Checkbox("Meshlet based Frustum Culling", &m_frustumCulling);
        ImGui::Checkbox("Screen Space Error based LoD selection", &m_screen_space_LoD);
        ImGui::Checkbox("Geo-Morphing", &m_geo_morphing);
        ImGui::Checkbox("Tessellation", &m_tessellation);
        if (m_screen_space_LoD)
            ImGui::InputFloat("Max error in pxl", &m_LoDScale, 0.01f, 1.0f, "%.2f");
        else 
            ImGui::InputFloat("LoD_0 Distance", &m_LoDScale, 0.01f, 1.0f, "%.2f");
        ImGui::SliderFloat("Debug Float", &m_debugFloatSlider, 0.0f, 10.0f, "%.2f");
        ImGui::SliderFloat("Displacement Scale", &m_displacementScale, 0.0f, 0.1f);
        ImGui::Checkbox("Triplanar Mapping instead of UV-coordiantes", &m_triplanarMapping);
        if (m_triplanarMapping)
        {
            ImGui::SliderFloat("Triplanar Texture Scale", &m_triplanarScale, 0.1f, 10.0f);
            ImGui::SliderFloat("Triplanar Blend Grade", &m_triplanarBlendGrade, 0.1f, 10.0f);
        }
        if (ImGui::Checkbox("Lock Camera Position Shader Constant", &m_lockCameraShaderConstant))
        {
            m_lockedCameraPos = m_cameraPos;
        }
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
            if (ImGui::Button("Reset Rotation"))
            {
                m_cameraPos = float3(-1.5, 0, 0);
                m_CameraRoll = 0.00;
                m_CameraYaw = PI / 2.0;
                m_autoRotationOffset = 0.0;
                //OutputDebugStringA(("X: " + std::to_string(m_cameraPos.x) + " Y: " + std::to_string(m_cameraPos.y) + " Z: " + std::to_string(m_cameraPos.z)+ " Roll: " + std::to_string(m_CameraRoll) + " Yaw: " + std::to_string(m_CameraYaw) + "\n").c_str());
            }
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
        ImGui::Checkbox("Normal based Lighting", &m_allowLighting);
        
        static int selected = static_cast<int>(m_shadingMode);  // Index of the selected option
        const char* options[] = { "Disable Debug Visals", "Show Meshlets", "Show Meshlet Grouping", "Show LoDs", "Show Tessellation Level", "Visualize Height Map"};
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

    // Performance Statistics window
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
        ImGui::Text("total  meshlet count: %d", (unsigned int)m_scene.m_total_meshlet_count);
        ImGui::Text("total   vertex count: %d", (unsigned int)m_scene.m_total_vertex_count);
        ImGui::Text("total triangle count: %d", (unsigned int)m_scene.m_totoal_triangle_count);
        ImGui::Separator();
        ImGui::Text("scene import   time: %.4f sec", m_scene.m_sceneImportTime);
        ImGui::Text("mesh parsing   time: %.4f sec", m_scene.m_modelParsingTime);
        ImGui::Text("hierachy-gen.  time: %.4f sec", m_scene.m_hierarchyGenTime);
        ImGui::Text("total pre-pro. time: %.4f sec", m_scene.m_totalPreProcessingTime);
        ImGui::Text("api setup      time: %.4f sec", m_apiSetupTime);
    }
    ImGui::End();

    // Help and Controls Window
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