#include <directx/d3dx12.h>
#include <MeshletLoD.h>

#include <psapi.h>

#include "HLSLnames.h"

#include "../shaders/Structures.fxh"

#include <Application.h>
#include <CommandQueue.h>
#include <Helpers.h>
#include <Window.h>

#include <imgui.h>
#include "imgui_impl_win32.h"
#include "imgui_impl_dx12.h"
#include "imgui_internal.h"

#include <wrl.h>
using namespace Microsoft::WRL;

#include <d3dcompiler.h>

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

MeshletLoD::MeshletLoD(const std::wstring& name, int width, int height, bool vSync)
    : super(name, width, height, vSync)
    , m_ScissorRect(CD3DX12_RECT(0, 0, LONG_MAX, LONG_MAX))
    , m_Viewport(CD3DX12_VIEWPORT(0.0f, 0.0f, static_cast<float>(width), static_cast<float>(height)))
    , m_FoV(90.0)
    , m_ContentLoaded(false)
{
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

void MeshletLoD::CreateBufferResource(
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

bool MeshletLoD::LoadContent()
{
    m_scene.init(m_model_file_path, 0);

    auto device = Application::Get().GetDevice();
    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_COPY);
    auto commandList = commandQueue->GetCommandList();



    D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
    heapDesc.NumDescriptors = (uint)m_scene.m_meshes.size() * 4 + 2;
    heapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    heapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    heapDesc.NodeMask = 0;
    device->CreateDescriptorHeap(&heapDesc, IID_PPV_ARGS(&m_CBV_SRV_UAV_Heap));

    unsigned int descriptorSize = device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV);
    CD3DX12_CPU_DESCRIPTOR_HANDLE srvHandle(m_CBV_SRV_UAV_Heap->GetCPUDescriptorHandleForHeapStart());
    
    // index buffers
    m_indexSrvHandle = m_CBV_SRV_UAV_Heap->GetGPUDescriptorHandleForHeapStart();

    std::vector<ComPtr<ID3D12Resource>> copyBuffers;
    for (int i = 0; i < m_scene.m_meshes.size(); i++)
    {
        // setup buffer
        ComPtr<ID3D12Resource> buffer;
        ComPtr<ID3D12Resource> copyBuffer;

        copyBuffers.push_back(copyBuffer);
        m_IndexBuffers.push_back(buffer);

        UpdateBufferResource(commandList, m_IndexBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)m_scene.m_meshes[i]->m_meshlet_vertices.size(), sizeof(uint), m_scene.m_meshes[i]->m_meshlet_vertices.data());

        // setup srv's
        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)m_scene.m_meshes[i]->m_meshlet_vertices.size();
        srvDesc.Buffer.StructureByteStride = sizeof(uint);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(m_IndexBuffers.back().Get(), &srvDesc, srvHandle);
        srvHandle.Offset(1, descriptorSize);
    }

    // vertex buffers
    m_vertexSrvHandle = m_indexSrvHandle;
    m_vertexSrvHandle.ptr += m_scene.m_meshes.size() * descriptorSize;

    for (int i = 0; i < m_scene.m_meshes.size(); i++)
    {
        // setup buffer
        ComPtr<ID3D12Resource> buffer;
        ComPtr<ID3D12Resource> copyBuffer;

        copyBuffers.push_back(copyBuffer);
        m_VertexBuffers.push_back(buffer);

        UpdateBufferResource(commandList, m_VertexBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)m_scene.m_meshes[i]->m_vertices.size(), sizeof(CustomVertex), m_scene.m_meshes[i]->m_vertices.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)m_scene.m_meshes[i]->m_vertices.size();
        srvDesc.Buffer.StructureByteStride = sizeof(CustomVertex);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(m_VertexBuffers.back().Get(), &srvDesc, srvHandle);
        srvHandle.Offset(1, descriptorSize);
    }

    // triangle buffers
    m_triangleSrvHandle = m_vertexSrvHandle;
    m_triangleSrvHandle.ptr += m_scene.m_meshes.size() * descriptorSize;


    for (int i = 0; i < m_scene.m_meshes.size(); i++)
    {
        // fill array with zeros to allign its size to 32bits
        while (((uint)m_scene.m_meshes[i]->m_meshlet_triangles.size()) % 4 != 0)
            m_scene.m_meshes[i]->m_meshlet_triangles.push_back(0);


        // setup buffer
        ComPtr<ID3D12Resource> buffer;
        ComPtr<ID3D12Resource> copyBuffer;

        copyBuffers.push_back(copyBuffer);
        m_TriangleBuffers.push_back(buffer);

        UpdateBufferResource(commandList, m_TriangleBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)m_scene.m_meshes[i]->m_meshlet_triangles.size() / 4, sizeof(uint), m_scene.m_meshes[i]->m_meshlet_triangles.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)m_scene.m_meshes[i]->m_meshlet_triangles.size() / 4; // store 4 chars per 32bit block
        srvDesc.Buffer.StructureByteStride = sizeof(uint);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(m_TriangleBuffers.back().Get(), &srvDesc, srvHandle);
        srvHandle.Offset(1, descriptorSize);
    }

    // draw task buffers
    m_drawTasksSrvHandle = m_triangleSrvHandle;
    m_drawTasksSrvHandle.ptr += m_scene.m_meshes.size() * descriptorSize;

    for (int i = 0; i < m_scene.m_meshes.size(); i++)
    {
        // setup buffer
        ComPtr<ID3D12Resource> buffer;
        ComPtr<ID3D12Resource> copyBuffer;

        copyBuffers.push_back(copyBuffer);
        m_DrawTasksBuffers.push_back(buffer);

        UpdateBufferResource(commandList, m_DrawTasksBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)m_scene.m_meshes[i]->m_draw_tasks.size(), sizeof(DrawTask), m_scene.m_meshes[i]->m_draw_tasks.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)m_scene.m_meshes[i]->m_draw_tasks.size();
        srvDesc.Buffer.StructureByteStride = sizeof(DrawTask);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(m_DrawTasksBuffers.back().Get(), &srvDesc, srvHandle);
        srvHandle.Offset(1, descriptorSize);
    }

    // objects buffer
    m_objectsSrvHandle = m_drawTasksSrvHandle;
    m_objectsSrvHandle.ptr += m_scene.m_meshes.size() * descriptorSize;

    ComPtr<ID3D12Resource> objectsCopyBuffer;
    UpdateBufferResource(commandList, m_ObjectsBuffer.GetAddressOf(), objectsCopyBuffer.GetAddressOf(),
        (uint)m_scene.m_scene_objects.size(), sizeof(SceneObject), m_scene.m_scene_objects.data());

    D3D12_SHADER_RESOURCE_VIEW_DESC objectsSrvDesc = {};
    objectsSrvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
    objectsSrvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    objectsSrvDesc.Buffer.FirstElement = 0;
    objectsSrvDesc.Buffer.NumElements = (uint)m_scene.m_scene_objects.size();
    objectsSrvDesc.Buffer.StructureByteStride = sizeof(SceneObject);
    objectsSrvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

    device->CreateShaderResourceView(m_ObjectsBuffer.Get(), &objectsSrvDesc, srvHandle);
    srvHandle.Offset(1, descriptorSize);

    // meshlet counts buffer
    m_meshletCountsSrvHandle = m_objectsSrvHandle;
    m_meshletCountsSrvHandle.ptr += descriptorSize;

    ComPtr<ID3D12Resource> meshletCountsCopyBuffer;
    UpdateBufferResource(commandList, m_MeshletCountsBuffer.GetAddressOf(), meshletCountsCopyBuffer.GetAddressOf(),
        (uint)m_scene.m_meshlet_counts.size(), sizeof(uint), m_scene.m_meshlet_counts.data());

    D3D12_SHADER_RESOURCE_VIEW_DESC meshletCountsSrvDesc = {};
    meshletCountsSrvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
    meshletCountsSrvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    meshletCountsSrvDesc.Buffer.FirstElement = 0;
    meshletCountsSrvDesc.Buffer.NumElements = (uint)m_scene.m_meshlet_counts.size();
    meshletCountsSrvDesc.Buffer.StructureByteStride = sizeof(uint);
    meshletCountsSrvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

    device->CreateShaderResourceView(m_MeshletCountsBuffer.Get(), &meshletCountsSrvDesc, srvHandle);
    srvHandle.Offset(1, descriptorSize);

    // Create the descriptor heap for the depth-stencil view.
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&m_DSVHeap)));

    // Load the pixel shader.
    ComPtr<ID3DBlob> pixelShaderBlob;
    ThrowIfFailed(D3DReadFileToBlob(L"PixelShader.cso", &pixelShaderBlob));

    // Load the mesh shader.
    ComPtr<ID3DBlob> meshShaderBlob;
    ThrowIfFailed(D3DReadFileToBlob(L"MeshShader.cso", &meshShaderBlob));

    // Load the task shader.
    ComPtr<ID3DBlob> taskShaderBlob;
    ThrowIfFailed(D3DReadFileToBlob(L"TaskShader.cso", &taskShaderBlob));


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

    // A single 32-bit constant root parameter that is used by the vertex shader.
    CD3DX12_ROOT_PARAMETER1 rootParameters[8];
    rootParameters[0].InitAsConstants(sizeof(Constants) / 4, 0, 0, D3D12_SHADER_VISIBILITY_ALL);

    // vertices
    CD3DX12_DESCRIPTOR_RANGE1 srvVertexRange;
    srvVertexRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_meshes.size(), 0, 1);
    rootParameters[2].InitAsDescriptorTable(1, &srvVertexRange, D3D12_SHADER_VISIBILITY_MESH);

    // indices
    CD3DX12_DESCRIPTOR_RANGE1 srvIndexRange;
    srvIndexRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_meshes.size(), 0, 2);
    rootParameters[1].InitAsDescriptorTable(1, &srvIndexRange, D3D12_SHADER_VISIBILITY_MESH);

    // triangles
    CD3DX12_DESCRIPTOR_RANGE1 srvTrianglesRange;
    srvTrianglesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_meshes.size(), 0, 3);
    rootParameters[4].InitAsDescriptorTable(1, &srvTrianglesRange, D3D12_SHADER_VISIBILITY_MESH);

    // draw tasks
    CD3DX12_DESCRIPTOR_RANGE1 srvDrawTasksRange;
    srvDrawTasksRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_meshes.size(), 0, 4);
    rootParameters[5].InitAsDescriptorTable(1, &srvDrawTasksRange, D3D12_SHADER_VISIBILITY_AMPLIFICATION);

    // objects
    rootParameters[6].InitAsShaderResourceView(0, 0);

    // meshlet counts
    rootParameters[7].InitAsShaderResourceView(1, 0);
    
    // indirect draw id
    rootParameters[3].InitAsConstants(1, 1, 0, D3D12_SHADER_VISIBILITY_ALL);
    

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


    // setup indirect draw attributes buffer
    ComPtr<ID3D12Resource> copy_indirectAtgumentBuffer;

    struct CommandStructure {
        unsigned int instanceID;
        D3D12_DISPATCH_MESH_ARGUMENTS dispatchArguments;
    };

    UpdateBufferResource(commandList,
        &m_indirectArgumentBuffer, &copy_indirectAtgumentBuffer,
        (uint)m_scene.m_indirect_attributes_withConstant.size() / 4, sizeof(CommandStructure), m_scene.m_indirect_attributes_withConstant.data());


    // transition to indirect argument state
    auto commandQueue2 = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_DIRECT);
    auto commandList2 = commandQueue2->GetCommandList();
    commandList2->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(m_indirectArgumentBuffer.Get(),
        D3D12_RESOURCE_STATE_COMMON, D3D12_RESOURCE_STATE_INDIRECT_ARGUMENT));


    // setup command signature
    D3D12_COMMAND_SIGNATURE_DESC commandSignatureDesc = {};
    D3D12_INDIRECT_ARGUMENT_DESC argumentDescs[2] = {};

    argumentDescs[0].Type = D3D12_INDIRECT_ARGUMENT_TYPE_CONSTANT;
    argumentDescs[0].Constant.RootParameterIndex = 3;
    argumentDescs[0].Constant.DestOffsetIn32BitValues = 0;
    argumentDescs[0].Constant.Num32BitValuesToSet = 1;

    argumentDescs[1].Type = D3D12_INDIRECT_ARGUMENT_TYPE_DISPATCH_MESH;

    commandSignatureDesc.pArgumentDescs = argumentDescs;
    commandSignatureDesc.NumArgumentDescs = _countof(argumentDescs);
    commandSignatureDesc.ByteStride = sizeof(CommandStructure);
    device->CreateCommandSignature(&commandSignatureDesc, m_RootSignature.Get(), IID_PPV_ARGS(&m_commandSignature));


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
    rasterizerDesc.CullMode = D3D12_CULL_MODE_BACK;

    pipelineStateStream.pRootSignature = m_RootSignature.Get();
    pipelineStateStream.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
    pipelineStateStream.MS = CD3DX12_SHADER_BYTECODE(meshShaderBlob.Get());
    pipelineStateStream.AS = CD3DX12_SHADER_BYTECODE(taskShaderBlob.Get());
    pipelineStateStream.PS = CD3DX12_SHADER_BYTECODE(pixelShaderBlob.Get());
    pipelineStateStream.DSVFormat = DXGI_FORMAT_D32_FLOAT;
    pipelineStateStream.RTVFormats = rtvFormats;
    pipelineStateStream.RasterizerState = rasterizerDesc;
    
    D3D12_PIPELINE_STATE_STREAM_DESC pipelineStateStreamDesc = {
        sizeof(PipelineStateStream), &pipelineStateStream
    };
    ThrowIfFailed(device->CreatePipelineState(&pipelineStateStreamDesc, IID_PPV_ARGS(&m_PipelineState)));

    auto fenceValue = commandQueue->ExecuteCommandList(commandList);
    commandQueue->WaitForFenceValue(fenceValue);

    fenceValue = commandQueue2->ExecuteCommandList(commandList2);
    commandQueue2->WaitForFenceValue(fenceValue);

    m_ContentLoaded = true;

    initImGui();

    // Resize/Create the depth buffer.
    ResizeDepthBuffer(GetClientWidth(), GetClientHeight());

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
}

void MeshletLoD::OnUpdate(UpdateEventArgs& e)
{
    static uint64_t frameCount = 0;
    static double totalTime = 0.0;

   
    m_frameTime = e.ElapsedTime;

    super::OnUpdate(e);

    totalTime += e.ElapsedTime;
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


    // Update the view matrix.
    const XMVECTOR eyePosition = XMVectorSet(0, 0, 0, 1);
    const XMVECTOR focusPoint = XMVectorSet(1, 0, 0, 1);
    const XMVECTOR upDirection = XMVectorSet(0, 1, 0, 0);
    m_ViewMatrix = XMMatrixLookAtLH(eyePosition, focusPoint, upDirection);
    
    // Camera Controll
    if (!m_freeCamera)
    {
        m_autoCameraDistance = std::clamp(m_autoCameraDistance - ImGui::GetIO().MouseWheel, 5.0f, 200.0f);
    }
    else
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

            m_CameraRoll = lerp(m_CameraRoll, (float)-PI * 0.1f, 2.5f * static_cast<float>(e.ElapsedTime));
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
    m_ProjectionMatrix = XMMatrixPerspectiveFovLH(XMConvertToRadians(m_FoV), aspectRatio, 0.1f, 1000.0f);

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

// Clear a render target.
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

    // Clear the render targets.
    {
        TransitionResource(commandList, backBuffer,
            D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET);

        ClearRTV(commandList, rtv, m_ClearColor);
        ClearDepth(commandList, dsv);
    }

    commandList->SetPipelineState(m_PipelineState.Get());
    commandList->SetGraphicsRootSignature(m_RootSignature.Get());

    commandList->RSSetViewports(1, &m_Viewport);
    commandList->RSSetScissorRects(1, &m_ScissorRect);

    commandList->OMSetRenderTargets(1, &rtv, FALSE, &dsv);

    // Update the MVP matrix
    /*
    XMMATRIX mvpMatrix = XMMatrixMultiply(m_ModelMatrix, m_ViewMatrix);
    mvpMatrix = XMMatrixMultiply(mvpMatrix, m_ProjectionMatrix);
    */
    XMMATRIX mvpMatrix = XMMatrixTranspose(XMMatrixMultiply(m_ViewMatrix, m_ProjectionMatrix));
    Constants constants;

    constants.ViewProjMat = mvpMatrix;
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
    constants.CurrTime = 0.0f;
    constants.BoolConstants = 0;
    if (m_frustumCulling) constants.BoolConstants |= FRUSTUM_CULLING_BIT_POS;
    if (m_coneCulling) constants.BoolConstants |= CONE_CULLING_BIT_POS;
    if (m_debugVisuals) constants.BoolConstants |= DEBUG_VISUALS_BIT_POS;



    commandList->SetGraphicsRoot32BitConstants(0, sizeof(Constants) / 4, &constants, 0);

    //commandList->DrawIndexedInstanced(_countof(g_Indicies), 1, 0, 0, 0);

    ID3D12DescriptorHeap* heaps[] = { m_CBV_SRV_UAV_Heap.Get()};
    commandList->SetDescriptorHeaps(_countof(heaps), heaps);

    // set bindless index buffers
    commandList->SetGraphicsRootDescriptorTable(1, m_indexSrvHandle);
    // set bindless vertex buffers
    commandList->SetGraphicsRootDescriptorTable(2, m_vertexSrvHandle);
    // set bindless triangles buffers
    commandList->SetGraphicsRootDescriptorTable(4, m_triangleSrvHandle);
    // set bindless draw tasks buffers
    commandList->SetGraphicsRootDescriptorTable(5, m_drawTasksSrvHandle);
    // set objects buffer
    commandList->SetGraphicsRootShaderResourceView(6, m_ObjectsBuffer.Get()->GetGPUVirtualAddress());
    // set meshlet counts buffer
    commandList->SetGraphicsRootShaderResourceView(7, m_MeshletCountsBuffer.Get()->GetGPUVirtualAddress());
    
    

    //commandList->DispatchMesh(1, 1, 1);

    commandList->ExecuteIndirect(
        m_commandSignature.Get(),
        (uint)m_scene.m_scene_objects.size(),
        m_indirectArgumentBuffer.Get(),
        0,
        nullptr,
        0
    );
    


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
    m_autoCameraDistance = std::clamp(m_autoCameraDistance - e.WheelDelta, 5.0f, 200.0f);
}

void MeshletLoD::initImGui() 
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();


    

    D3D12_DESCRIPTOR_HEAP_DESC descHeap = {};
    descHeap.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
    descHeap.NumDescriptors = 1;
    descHeap.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
    descHeap.NodeMask = 0;

    ThrowIfFailed(Application::Get().GetDevice().Get()->CreateDescriptorHeap(
        &descHeap, IID_PPV_ARGS(&m_ImGuiDescriptorHeap)
    ));

    
    // Setup Platform/Renderer backends
    ImGui_ImplWin32_Init(m_pWindow->GetWindowHandle());
    ImGui_ImplDX12_Init(Application::Get().GetDevice().Get(), 3,
        DXGI_FORMAT_R8G8B8A8_UNORM, m_ImGuiDescriptorHeap.Get(),
        m_ImGuiDescriptorHeap->GetCPUDescriptorHandleForHeapStart(),
        m_ImGuiDescriptorHeap->GetGPUDescriptorHandleForHeapStart()
    );   

}

void MeshletLoD::updateImGui()
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
    fpsHistoryOffset += fpsHistorySpeed * (float)m_frameTime;
    if (fpsHistoryOffset >= fpsHistorySize) fpsHistoryOffset = 0.0f;



    // Start the Dear ImGui frame
    ImGui_ImplDX12_NewFrame();
    ImGui_ImplWin32_NewFrame();
    ImGui::NewFrame();


    ImGui::SetNextWindowSize(ImVec2(450, 0));
    ImGui::SetNextWindowPos(ImVec2(12, 12), ImGuiCond_Always);
    ImGui::Begin("Settings");                           

    if (!ImGui::CollapsingHeader("Perfromance"))
    {
        ImGui::Text("%.1f FPS", m_fps);
        ImGui::Text("%.2fms Frame Time", m_frameTime * 1000);
        ImGui::SameLine();
        ImGui::PlotLines("\0", fpsHistory, fpsHistorySize);
        if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&memoryInfo, sizeof(memoryInfo))) 
            ImGui::Text("%.3fGB RAM Usage", (memoryInfo.WorkingSetSize) / 1073741824.0f);
        ImGui::Text("%.3fGB VRAM Usage", videoMemoryInfo.CurrentUsage / 1000000000.0f);   
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Render Settings"))
    {
        ImGui::Checkbox("Frustum Culling", &m_frustumCulling);
        ImGui::Checkbox("Cone Culling", &m_coneCulling);
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
            ImGui::SliderFloat("Camera speed", &m_cameraSpeed, 0.1f, 100.0f);
        }
        else
        {
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
         
        }
        static int lod = 0;
        ImGui::SliderInt("Level of Detail", &lod, 0, 10);
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Debug Settings"))
    {
        ImGui::Checkbox("Wireframe", &m_wireframe);
        ImGui::Checkbox("Debug Visuals", &m_debugVisuals);

        
        if (ImGui::Button("Force Lag"))
            Sleep(1000);
        
    }
    ImGui::Spacing();

    ImGui::End();


    ImGui::SetNextWindowSize(ImVec2(275, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 275 - 12), 12.0f), ImGuiCond_Always);
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