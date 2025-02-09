#include <directx/d3dx12.h>
#include <MeshletLoD.h>

#include <psapi.h>

#include "HLSLnames.h"

#include "../shaders/Structures.fxh"

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

void MeshletLoD::CreatePSO()
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

void MeshletLoD::CreateCullingPSO()
{
    m_ObjectCulling_PipelineState.Reset();
    auto device = Application::Get().GetDevice();

    struct PipelineStateStream
    {
        CD3DX12_PIPELINE_STATE_STREAM_ROOT_SIGNATURE pRootSignature;
        CD3DX12_PIPELINE_STATE_STREAM_CS CS;
    } pipelineStateStream;

    pipelineStateStream.pRootSignature = m_ObjectCulling_RootSignature.Get();
    pipelineStateStream.CS = CD3DX12_SHADER_BYTECODE(m_objectCullingComputeShaderBlob.Get());

    D3D12_PIPELINE_STATE_STREAM_DESC pipelineStateStreamDesc = {
        sizeof(PipelineStateStream), &pipelineStateStream
    };
    ThrowIfFailed(device->CreatePipelineState(&pipelineStateStreamDesc, IID_PPV_ARGS(&m_ObjectCulling_PipelineState)));
}


void recursiveTestHelper(aiNode* currentNode, float4x4 parentTransform, std::map<std::string, float4x4>* boneMap, uint recusrionCounter)
{
    std::string nodeName = currentNode->mName.data;
    float4x4 currentTransform = MeshletMesh::assimpToFloat4x4Matrix(currentNode->mTransformation);
    OutputDebugStringA((std::string("\nRecursion Depth: ") + std::to_string(recusrionCounter)).c_str());
    OutputDebugStringA((std::string("\nCurrent Node: ") + nodeName).c_str());
    OutputDebugStringA((std::string("\nLocal Transform:\n") + MatrixToString(currentTransform)).c_str());
    currentTransform = currentTransform * parentTransform;
    

    if (boneMap->find(nodeName) != boneMap->end())
    {
        OutputDebugStringA((std::string("Is a Bone!")).c_str());
        OutputDebugStringA((std::string("\nBone Space Transform (local * parent):\n") + MatrixToString(currentTransform)).c_str());
        float4x4 offset = (*boneMap)[nodeName];
        OutputDebugStringA((std::string("Inverse Bind Pose Transform:\n") + MatrixToString(offset)).c_str());
        OutputDebugStringA((std::string("Offset * globalTransform:\n") + MatrixToString(offset * currentTransform)).c_str());
    }

    recusrionCounter++;
    for (uint i = 0; i < currentNode->mNumChildren; i++)
    {
        recursiveTestHelper(currentNode->mChildren[i], currentTransform, boneMap, recusrionCounter);
    }
}


bool MeshletLoD::LoadContent()
{
    // transformation testing
    
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(m_model_file_path, aiProcess_Triangulate);// | aiProcess_MakeLeftHanded | aiProcess_FlipWindingOrder);
    scene->mMeshes[0]->mBones;

    scene->mRootNode;
    
    std::map<std::string, float4x4> boneMap;

    for (int b = 0; b < scene->mMeshes[0]->mNumBones; b++)
    {
        std::string boneName = scene->mMeshes[0]->mBones[b]->mName.C_Str();
        float4x4 transform_matrix = MeshletMesh::assimpToFloat4x4Matrix(scene->mMeshes[0]->mBones[b]->mOffsetMatrix);
        boneMap[boneName] = transform_matrix;
    }

    //OutputDebugStringA(("\n\nMesh NumBones: " + std::to_string(scene->mMeshes[0]->mNumBones) + " Animation NumChannels: " + std::to_string(scene->mAnimations[0]->mChannels[0]. ->mNumMeshChannels) + "\n\n").c_str());
    //OutputDebugStringA(("\n\nMesh Name: " + std::string(scene->mMeshes[0]->mName.C_Str()) + " Animation Channel 0 name: " + std::string(scene->mAnimations[0]->mMeshChannels[0]->mName.C_Str()) + "\n\n").c_str());

   
    //recursiveTestHelper(scene->mRootNode, XMMatrixIdentity(), &boneMap, 0);
    
    

    //OutputDebugStringA((std::string("Test Transform:\n") + MatrixToString(XMMatrixTranslation(10, 2, 0))).c_str());
    //OutputDebugStringA((std::string("Test Transform:\n") + MatrixToString(XMMatrixRotationQuaternion(XMVectorSet(0, 0, 0, 1)))).c_str());
    //OutputDebugStringA((std::string("Test Transform2:\n") + MatrixToString(XMMatrixRotationQuaternion(XMVectorSet(0.707107, 0, 0, 0.707107)))).c_str());
    //OutputDebugStringA((std::string("Test Transform3:\n") + MatrixToString(XMMatrixRotationQuaternion(XMVectorSet(-0.325058, -0.000000, -0.325058, 0.888074)))).c_str());
    
    /*
    float3 p1(0, 0, 0);
    float3 p2(1, 1, 1);
    float scaleFactor = -0.2;

    XMVECTOR vecTest = XMVectorLerp(DirectX::XMLoadFloat3(&(p1)), DirectX::XMLoadFloat3(&(p2)), scaleFactor);
    OutputDebugStringA(("\nTest VEC LERP" + MatrixToString(XMMatrixTranslationFromVector(vecTest))).c_str());
    */

    m_scene.init(m_model_file_path, 0);
    
 

    



    /*
    OutputDebugStringA((std::string("\nBone Location Keyframes:\n")).c_str());

    for (auto bone : m_scene.animator.m_CurrentAnimation->m_Bones)
    {
        OutputDebugStringA(("\nBone Name: \"" + bone.m_Name + "\" contains " + std::to_string(bone.m_Positions.size()) + " position keyframes").c_str());
        for (auto key : bone.m_Positions)
        {
            OutputDebugStringA(("\nAnimation time: " + std::to_string(key.timeStamp)).c_str());
            OutputDebugStringA(("\nOffset: (" + std::to_string(key.position.x) + ", " + std::to_string(key.position.y) + ", " + std::to_string(key.position.x) + ")\n").c_str());
        }
        OutputDebugStringA(("\nBone Name: \"" + bone.m_Name + "\" contains " + std::to_string(bone.m_Positions.size()) + " rotation keyframes").c_str());
        for (auto key : bone.m_Rotations)
        {
            OutputDebugStringA(("\nAnimation time: " + std::to_string(key.timeStamp)).c_str());
            OutputDebugStringA(("\nOrientation: (" + std::to_string(XMVectorGetX(key.orientation)) + ", " + std::to_string(XMVectorGetY(key.orientation)) + ", " + std::to_string(XMVectorGetZ(key.orientation)) + ", " + std::to_string(XMVectorGetW(key.orientation)) + ")\n").c_str());
        }
    }
    */

    auto device = Application::Get().GetDevice();
    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_COPY);
    auto commandList = commandQueue->GetCommandList();



    D3D12_DESCRIPTOR_HEAP_DESC heapDesc = {};
    heapDesc.NumDescriptors = (uint)m_scene.m_meshes.size() * 4 + 4 + (uint)m_scene.m_preBakedAnimations.size(); // 4 bindless buffers that scale with mesh count + 1 bindless with animations + 4 single buffers
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

    // visible object count buffer
    m_visibleObjectCountSrvHandle = m_meshletCountsSrvHandle;
    m_visibleObjectCountSrvHandle.ptr += descriptorSize;
    uint initialObjectCount = (uint)m_scene.m_scene_objects.size();

    ComPtr<ID3D12Resource> visibleObjectCountCopyBuffer;
    UpdateBufferResource(commandList, m_objectCountBuffer.GetAddressOf(), visibleObjectCountCopyBuffer.GetAddressOf(),
        1, sizeof(uint), &initialObjectCount);

    D3D12_SHADER_RESOURCE_VIEW_DESC visibleObjectCountSrvDesc = {};
    visibleObjectCountSrvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
    visibleObjectCountSrvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    visibleObjectCountSrvDesc.Buffer.FirstElement = 0;
    visibleObjectCountSrvDesc.Buffer.NumElements = 1;
    visibleObjectCountSrvDesc.Buffer.StructureByteStride = sizeof(uint);
    visibleObjectCountSrvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

    device->CreateShaderResourceView(m_MeshletCountsBuffer.Get(), &visibleObjectCountSrvDesc, srvHandle);
    srvHandle.Offset(1, descriptorSize);


    // bone transformation matrices buffers
    m_boneMatricesSrvHandle = m_visibleObjectCountSrvHandle;
    m_boneMatricesSrvHandle.ptr += descriptorSize;

    for (int i = 0; i < m_scene.m_preBakedAnimations.size(); i++)
    {
        // setup buffer
        ComPtr<ID3D12Resource> buffer;
        ComPtr<ID3D12Resource> copyBuffer;

        copyBuffers.push_back(copyBuffer);
        m_BoneMatricesBuffers.push_back(buffer);

        UpdateBufferResource(commandList, m_BoneMatricesBuffers.back().GetAddressOf(), copyBuffers.back().GetAddressOf(),
            (uint)m_scene.m_preBakedAnimations[i]->matrices.size(), sizeof(float4x4), m_scene.m_preBakedAnimations[i]->matrices.data());

        D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
        srvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
        srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
        srvDesc.Buffer.FirstElement = 0;
        srvDesc.Buffer.NumElements = (uint)m_scene.m_preBakedAnimations[i]->matrices.size();
        srvDesc.Buffer.StructureByteStride = sizeof(float4x4);
        srvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

        device->CreateShaderResourceView(m_BoneMatricesBuffers.back().Get(), &srvDesc, srvHandle);
        srvHandle.Offset(1, descriptorSize);
    }

    // animations meta data buffer
    m_animationMetaDataSrvHandle = m_boneMatricesSrvHandle;
    m_animationMetaDataSrvHandle.ptr += m_scene.m_preBakedAnimations.size() * descriptorSize;

    assert((uint)m_scene.m_animationMetaData.size() == (uint)m_scene.m_preBakedAnimations.size());

    ComPtr<ID3D12Resource> animationMetaDataCopyBuffer;
    if ((uint)m_scene.m_animationMetaData.size())
        UpdateBufferResource(commandList, m_AnimationMetaDataBuffer.GetAddressOf(), animationMetaDataCopyBuffer.GetAddressOf(),
            (uint)m_scene.m_animationMetaData.size(), sizeof(AnimationMetaData), m_scene.m_animationMetaData.data());

    D3D12_SHADER_RESOURCE_VIEW_DESC animationMetaDataSrvDesc = {};
    animationMetaDataSrvDesc.ViewDimension = D3D12_SRV_DIMENSION_BUFFER;
    animationMetaDataSrvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
    animationMetaDataSrvDesc.Buffer.FirstElement = 0;
    animationMetaDataSrvDesc.Buffer.NumElements = (uint)m_scene.m_animationMetaData.size();
    animationMetaDataSrvDesc.Buffer.StructureByteStride = sizeof(AnimationMetaData);
    animationMetaDataSrvDesc.Buffer.Flags = D3D12_BUFFER_SRV_FLAG_NONE;

    device->CreateShaderResourceView(m_AnimationMetaDataBuffer.Get(), &animationMetaDataSrvDesc, srvHandle);
    srvHandle.Offset(1, descriptorSize);


    // Create the descriptor heap for the depth-stencil view.
    D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
    dsvHeapDesc.NumDescriptors = 1;
    dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
    dsvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_NONE;
    ThrowIfFailed(device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&m_DSVHeap)));

    // Load the pixel shader.
    
    ThrowIfFailed(D3DReadFileToBlob(L"PixelShader.cso", &m_pixelShaderBlob));

    // Load the mesh shader.
    ThrowIfFailed(D3DReadFileToBlob(L"MeshShader.cso", &m_meshShaderBlob));

    // Load the task shader.
    ThrowIfFailed(D3DReadFileToBlob(L"TaskShader.cso", &m_taskShaderBlob));

    // Load the object culling compute shader.
    ThrowIfFailed(D3DReadFileToBlob(L"ObjectCulling_ComputeShader.cso", &m_objectCullingComputeShaderBlob));


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


    CD3DX12_ROOT_PARAMETER1 rootParameters[10];
    // constants
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

    // skeletal animation bone matrices
    CD3DX12_DESCRIPTOR_RANGE1 srvBoneMatricesRange;
    srvBoneMatricesRange.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, (uint)m_scene.m_meshes.size(), 0, 5);
    rootParameters[8].InitAsDescriptorTable(1, &srvBoneMatricesRange, D3D12_SHADER_VISIBILITY_MESH);

    // animation meta data
    rootParameters[9].InitAsShaderResourceView(2, 0);
    

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




    // setup object culling root signature
    CD3DX12_ROOT_PARAMETER1 objectCullingRootParameters[6];

    // constants
    objectCullingRootParameters[0].InitAsConstants(sizeof(Constants) / 4, 0, 0, D3D12_SHADER_VISIBILITY_ALL);
    objectCullingRootParameters[5].InitAsConstants(sizeof(uint) / 4, 1, 0, D3D12_SHADER_VISIBILITY_ALL);

    // objects
    objectCullingRootParameters[1].InitAsShaderResourceView(0, 0);

    // meshlet counts
    objectCullingRootParameters[2].InitAsShaderResourceView(1, 0);

    // visible object count buffer
    objectCullingRootParameters[3].InitAsUnorderedAccessView(0, 0);

    // inirect arguments buffer
    objectCullingRootParameters[4].InitAsUnorderedAccessView(1, 0);

    CD3DX12_VERSIONED_ROOT_SIGNATURE_DESC objectCullingRootSignatureDescription;
    objectCullingRootSignatureDescription.Init_1_1(_countof(objectCullingRootParameters), objectCullingRootParameters, 0, nullptr, D3D12_ROOT_SIGNATURE_FLAG_NONE);

    // Serialize the root signature.
    ComPtr<ID3DBlob> objectCullingRootSignatureBlob;
    ComPtr<ID3DBlob> objectCullingErrorBlob;
    ThrowIfFailed(D3DX12SerializeVersionedRootSignature(&objectCullingRootSignatureDescription,
        featureData.HighestVersion, &objectCullingRootSignatureBlob, &objectCullingErrorBlob));

    // Create the root signature.
    ThrowIfFailed(device->CreateRootSignature(0, objectCullingRootSignatureBlob->GetBufferPointer(),
        objectCullingRootSignatureBlob->GetBufferSize(), IID_PPV_ARGS(&m_ObjectCulling_RootSignature)));




    // setup indirect draw attributes buffer
    ComPtr<ID3D12Resource> copy_indirectAtgumentBuffer;

    UpdateBufferResource(commandList,
        &m_indirectArgumentBuffer, &copy_indirectAtgumentBuffer,
        (uint)m_scene.m_indirect_attributes.size(), sizeof(CommandStructure), m_scene.m_indirect_attributes.data());


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

    CreatePSO();
    CreateCullingPSO();

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

    m_RootSignature.Reset();
    m_PipelineState.Reset();

    m_ObjectCulling_RootSignature.Reset();
    m_ObjectCulling_PipelineState.Reset();

    // clear & reset model buffers
    m_CBV_SRV_UAV_Heap.Reset();

    m_IndexBuffers.clear();
    m_VertexBuffers.clear();
    m_TriangleBuffers.clear();
    m_DrawTasksBuffers.clear();
    m_ObjectsBuffer.Reset();
    m_MeshletCountsBuffer.Reset();
    m_objectCountBuffer.Reset();
    m_indirectArgumentBuffer.Reset();

    m_commandSignature.Reset();
}

void MeshletLoD::OnUpdate(UpdateEventArgs& e)
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
    m_ProjectionMatrix = XMMatrixPerspectiveFovLH(XMConvertToRadians(m_FoV), aspectRatio, 0.1f, 100000.0f);

    updateImGui();


    
}

void MeshletLoD::UpdateObjectsBuffer()
{
    Application::Get().Flush();
    auto device = Application::Get().GetDevice();
    auto commandQueue = Application::Get().GetCommandQueue(D3D12_COMMAND_LIST_TYPE_COPY);
    auto commandList = commandQueue->GetCommandList();
    

    
    
    ComPtr<ID3D12Resource> copyBuffer;

    

    //UpdateBufferResource(commandList, m_BoneMatricesBuffers[i].GetAddressOf(), copyBuffers.back().GetAddressOf(),
    //    (uint)m_scene.animator.m_FinalBoneMatrices.size(), sizeof(float4x4), m_scene.animator.m_FinalBoneMatrices.data());//(uint)m_scene.m_meshes[i]->m_draw_tasks.size(), sizeof(DrawTask), m_scene.m_meshes[i]->m_draw_tasks.data());

    size_t bufferSize = (uint)m_scene.m_scene_objects.size() * sizeof(SceneObject);

    ThrowIfFailed(device->CreateCommittedResource(
        &CD3DX12_HEAP_PROPERTIES(D3D12_HEAP_TYPE_UPLOAD),
        D3D12_HEAP_FLAG_NONE,
        &CD3DX12_RESOURCE_DESC::Buffer(bufferSize),
        D3D12_RESOURCE_STATE_GENERIC_READ,
        nullptr,
        IID_PPV_ARGS(copyBuffer.GetAddressOf())));

    D3D12_SUBRESOURCE_DATA subresourceData = {};
    subresourceData.pData = m_scene.m_scene_objects.data();
    subresourceData.RowPitch = bufferSize;
    subresourceData.SlicePitch = subresourceData.RowPitch;

    UpdateSubresources(commandList.Get(),
        *m_ObjectsBuffer.GetAddressOf(), *copyBuffer.GetAddressOf(),
        0, 0, 1, &subresourceData);
    

    auto fenceValue = commandQueue->ExecuteCommandList(commandList);
    commandQueue->WaitForFenceValue(fenceValue);
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


    // Update the MVP matrix
    XMMATRIX mvpMatrix = XMMatrixTranspose(m_ViewMatrix * m_ProjectionMatrix);
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
    constants.CurrTime = m_totalRunTime;
    constants.BoolConstants = 0;
    if (m_frustumCulling) constants.BoolConstants |= FRUSTUM_CULLING_BIT_POS;
    if (m_coneCulling) constants.BoolConstants |= CONE_CULLING_BIT_POS;
    if (m_debugVisuals) constants.BoolConstants |= DEBUG_VISUALS_BIT_POS;


    // Clear the render targets.
    {
        TransitionResource(commandList, backBuffer,
            D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET);

        ClearRTV(commandList, rtv, m_ClearColor);
        ClearDepth(commandList, dsv);
    }




    // object culling compute pass
    commandList->SetPipelineState(m_ObjectCulling_PipelineState.Get());
    commandList->SetComputeRootSignature(m_ObjectCulling_RootSignature.Get());

    commandList->SetComputeRoot32BitConstants(0, sizeof(Constants) / 4, &constants, 0);
    uint maxObjectCount = (uint)m_scene.m_scene_objects.size();
    commandList->SetComputeRoot32BitConstants(5, sizeof(uint) / 4, &maxObjectCount, 0);

    ID3D12DescriptorHeap* heaps[] = { m_CBV_SRV_UAV_Heap.Get() };
    commandList->SetDescriptorHeaps(_countof(heaps), heaps);

    // set objects buffer
    commandList->SetComputeRootShaderResourceView(1, m_ObjectsBuffer.Get()->GetGPUVirtualAddress());
    // set meshlet counts buffer
    commandList->SetComputeRootShaderResourceView(2, m_MeshletCountsBuffer.Get()->GetGPUVirtualAddress());
    // set visible object count buffer
    commandList->SetComputeRootUnorderedAccessView(3, m_objectCountBuffer.Get()->GetGPUVirtualAddress());
    // set indirect arguments buffer
    commandList->SetComputeRootUnorderedAccessView(4, m_indirectArgumentBuffer.Get()->GetGPUVirtualAddress());
    

    uint32_t groupCountX = (maxObjectCount + GROUP_SIZE - 1) / GROUP_SIZE;
    if (m_objectCulling)
        commandList->Dispatch(groupCountX, 1, 1);



    // main render pass
    commandList->SetPipelineState(m_PipelineState.Get());
    commandList->SetGraphicsRootSignature(m_RootSignature.Get());

    commandList->RSSetViewports(1, &m_Viewport);
    commandList->RSSetScissorRects(1, &m_ScissorRect);

    commandList->OMSetRenderTargets(1, &rtv, FALSE, &dsv);

    commandList->SetGraphicsRoot32BitConstants(0, sizeof(Constants) / 4, &constants, 0);

    //commandList->DrawIndexedInstanced(_countof(g_Indicies), 1, 0, 0, 0);

    /*
    ID3D12DescriptorHeap* heaps[] = { m_CBV_SRV_UAV_Heap.Get()};
    commandList->SetDescriptorHeaps(_countof(heaps), heaps);
    */

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
    // set bone matrices buffers
    commandList->SetGraphicsRootDescriptorTable(8, m_boneMatricesSrvHandle);
    // set animation meta data buffer
    if ((uint)m_scene.m_animationMetaData.size())
        commandList->SetGraphicsRootShaderResourceView(9, m_AnimationMetaDataBuffer.Get()->GetGPUVirtualAddress());
    
    

    //commandList->DispatchMesh(1, 1, 1);

    commandList->ExecuteIndirect(
        m_commandSignature.Get(),
        (uint)m_scene.m_scene_objects.size(),
        m_indirectArgumentBuffer.Get(),
        0,
        m_objectCountBuffer.Get(),
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
    m_autoCameraDistance = std::clamp(m_autoCameraDistance - e.WheelDelta, 1.0f, 200.0f);
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
    fpsHistoryOffset += clamp(fpsHistorySpeed * (float)m_frameTime, 0.0f, 1.0f);
    if (fpsHistoryOffset >= fpsHistorySize) fpsHistoryOffset = 0.0f;



    // Start the Dear ImGui frame
    
    ImGui_ImplWin32_NewFrame();
    ImGui_ImplDX12_NewFrame();
    ImGui::NewFrame();


    ImGui::SetNextWindowSize(ImVec2(375, 0));
    ImGui::SetNextWindowPos(ImVec2(275, 12), ImGuiCond_Always);
    ImGui::Begin("Settings");                           

    

    if (!ImGui::CollapsingHeader("Render Settings"))
    {
        ImGui::Checkbox("Meshlet based Frustum Culling", &m_frustumCulling);
        ImGui::Checkbox("Meshlet based Cone Culling", &m_coneCulling);
        ImGui::Checkbox("Object Culling", &m_objectCulling);
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
            ImGui::SliderFloat("Camera speed", &m_cameraSpeed, 0.1f, 1000.0f);
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
            Application::Get().Flush();
            UnloadContent();
            LoadContent();
        }
        static int lod = 0;
        ImGui::SliderInt("Level of Detail", &lod, 0, 10);
    }
    ImGui::Spacing();

    if (!ImGui::CollapsingHeader("Debug Settings"))
    {
        if (ImGui::Checkbox("Wireframe", &m_wireframe)) CreatePSO();
        ImGui::SameLine();
        if (ImGui::Checkbox("Back Face Culling", &m_backFaceCulling)) CreatePSO();
        ImGui::SameLine();
        if (ImGui::Checkbox("Debug Visuals", &m_debugVisuals)) CreatePSO();
        
        if (ImGui::Button("Force Lag"))
            Sleep(1000);
        
    }
    ImGui::Spacing();

    ImGui::End();

    // Animation Window
    ImGui::SetNextWindowSize(ImVec2(275, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 275 - 12), 12.0f), ImGuiCond_Always);
    ImGui::Begin("Animation");

    if (!ImGui::CollapsingHeader("General Animation Settings"))
    {
        if (ImGui::Button("Reset All"));
        if (ImGui::Button("Randomize Offsets"));
    }

    ImGui::Dummy(ImVec2(0.0f, 7.0f));

    if (!ImGui::CollapsingHeader("Individual Object Animation"))
    {
        static int currentObjectSelection = 0; // Currently selected object index
        ImGui::Text("Current Object Selection:");
        ImGui::SetNextItemWidth(260);
        ImGui::Combo("", &currentObjectSelection, m_scene.m_sceneObjectNamesCharP.data(), m_scene.m_sceneObjectNamesCharP.size());

        int currentAnimationSelection = m_scene.m_scene_objects[currentObjectSelection].animation_id;
        for (uint a = 0; a < m_scene.m_meshes[m_scene.m_scene_objects[currentObjectSelection].mesh_id]->m_animations.size(); a++)
        {
            if (m_scene.m_scene_objects[currentObjectSelection].animation_id == m_scene.m_meshes[m_scene.m_scene_objects[currentObjectSelection].mesh_id]->m_animations[a].totalAnimationIndex)
                currentAnimationSelection = a;
        }
        
        currentAnimationSelection++;

        ImGui::Text("Set Animation:");
        ImGui::SetNextItemWidth(260);
        ImGui::Combo("##xx", &currentAnimationSelection,
            m_scene.m_meshes[m_scene.m_scene_objects[currentObjectSelection].mesh_id]->m_animationNamesCharP.data(),
            m_scene.m_meshes[m_scene.m_scene_objects[currentObjectSelection].mesh_id]->m_animationNamesCharP.size());

        if (currentAnimationSelection == 0)
            m_scene.m_scene_objects[currentObjectSelection].animation_id = -1;
        else
            m_scene.m_scene_objects[currentObjectSelection].animation_id = m_scene.m_meshes[m_scene.m_scene_objects[currentObjectSelection].mesh_id]->m_animations[currentAnimationSelection - 1].totalAnimationIndex;

        if (ImGui::Button("Disable Animation"))
        {
            m_scene.m_scene_objects[currentObjectSelection].animation_id = -1;
        }

        ImGui::Text("Animation Time Offset:");
        ImGui::SetNextItemWidth(160);
        ImGui::InputFloat("(ctrl) - 0.01s", &m_scene.m_scene_objects[currentObjectSelection].animation_time_offset, 1.0f, 0.01f, "%.2f");

        ImGui::Text("Animation Speed:");
        ImGui::SetNextItemWidth(160);
        ImGui::InputFloat("(ctrl) - 1%", &m_scene.m_scene_objects[currentObjectSelection].animation_speed, 1.0f, 0.01f, "%.2f");
    }
    
    ImGui::Dummy(ImVec2(0.0f, 7.0f));
    ImGui::Separator();
    ImGui::Dummy(ImVec2(0.0f, 3.0f));
    ImGui::SetCursorPosX(ImGui::GetContentRegionMax().x - ImGui::CalcTextSize("Apply Changes").x - ImGui::GetStyle().FramePadding.x * 2);
    if (ImGui::Button("Apply Changes")) 
    {
        UpdateObjectsBuffer();
    }

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
        ImGui::Text("unique  meshes count: %d", (unsigned int)m_scene.m_meshes.size());
        ImGui::Text("total  meshlet count: %d", (unsigned int)m_scene.m_draw_task_count);
        ImGui::Text("total   vertex count: %d", (unsigned int)m_scene.m_vertex_count);
        ImGui::Text("total triangle count: %d", (unsigned int)m_scene.m_triangles_count);

        if (!ImGui::CollapsingHeader("Scene Processing Times"))
        {
            ImGui::Text("model load time:         %.4f sec", m_scene.m_modelLoadTime.count());
            ImGui::Text("LoD generation time:     %.4f sec", m_scene.m_totalLoDGenTime.count());
            ImGui::Text("meshlet generation time: %.4f sec", m_scene.m_totalMeshletGenTime.count());
            ImGui::Separator();
            ImGui::Text("total load time:         %.4f sec", m_scene.m_totalTime.count());
        }

    }

    ImGui::End();


    // Help Window
    ImGui::SetNextWindowSize(ImVec2(275, 0));
    ImGui::SetNextWindowPos(ImVec2((float)(GetClientWidth() - 275 - 275 - 12 - 3), 12.0f), ImGuiCond_Always);
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