#pragma once

#include <Game.h>
#include <Window.h>
#include <vector>

#include <DirectXMath.h>

#include "Scene.h"

class MeshletLoD : public Game
{
public:
    using super = Game;

    // loads and initializes resources needed for rendering
    MeshletLoD(const std::wstring& name, int width, int height, bool vSync = false);
    virtual bool LoadContent() override;

    // frees resources allocted by the class
    virtual void UnloadContent() override;

protected:

    // updates camera and UI aswell as their related render variables
    // will run every frame and before the OnRender() function
    virtual void OnUpdate(UpdateEventArgs& e) override;

    // updates constants buffer and renders the current frame
    virtual void OnRender(RenderEventArgs& e) override;

    // handles keyboard input by user
    virtual void OnKeyPressed(KeyEventArgs& e) override;

    // handles mouse wheel input by user
    virtual void OnMouseWheel(MouseWheelEventArgs& e) override;

    // handles rescaling of apllication window like updating render targets
    virtual void OnResize(ResizeEventArgs& e) override; 

private:

    // ImGui related functions
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_ImGuiDescriptorHeap;
    void initImGui();
    void updateImGui();

    // application load helper functions
    void createPSO();
    
    void setupConstantsUploadBuffer();
    
    template <typename T>
    void setupSrvAndBuffer(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                           D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle,
                           D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                           std::vector<T>& cpuBuffer,
                           Microsoft::WRL::ComPtr<ID3D12Resource>& gpuBuffer,
                           std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                           Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                           unsigned int descriptorSize, D3D12_RESOURCE_FLAGS bufferFlags);
   
    template <typename T>
    void setupBindlessSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                                    D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle,
                                    D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                    std::vector<std::vector<T>>& cpuBuffers,
                                    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                    Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                    unsigned int descriptorSize);
    
    void setupBindlessUCharToUIntSrvAndBuffers(D3D12_GPU_DESCRIPTOR_HANDLE& srvGpuHandle,
                                               D3D12_GPU_DESCRIPTOR_HANDLE& nextAvailableGpuSrvHandle,
                                               D3D12_CPU_DESCRIPTOR_HANDLE& nextAvailableCpuSrvHandle,
                                               std::vector<std::vector<unsigned char>>& unsignedCharsCpuBuffers,
                                               std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& gpuBuffers,
                                               std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>>& copyBuffers,
                                               Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7>& commandList,
                                               unsigned int descriptorSize);

    void TransitionResource(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
        Microsoft::WRL::ComPtr<ID3D12Resource> resource,
        D3D12_RESOURCE_STATES beforeState, D3D12_RESOURCE_STATES afterState);

    // Clear a render target view.
    void ClearRTV(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
        D3D12_CPU_DESCRIPTOR_HANDLE rtv, FLOAT* clearColor);

    // Clear the depth of a depth-stencil view.
    void ClearDepth(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
        D3D12_CPU_DESCRIPTOR_HANDLE dsv, FLOAT depth = 1.0f );

    // Create a GPU buffer.
    void UpdateBufferResource(Microsoft::WRL::ComPtr<ID3D12GraphicsCommandList7> commandList,
        ID3D12Resource** pDestinationResource, ID3D12Resource** pIntermediateResource,
        size_t numElements, size_t elementSize, const void* bufferData, 
        D3D12_RESOURCE_FLAGS flags = D3D12_RESOURCE_FLAG_NONE );

    // Resize the depth buffer to match the size of the client area.
    void ResizeDepthBuffer(int width, int height);
    
    uint64_t m_FenceValues[Window::BufferCount] = {};

    // Depth buffer.
    Microsoft::WRL::ComPtr<ID3D12Resource> m_DepthBuffer;
    // Descriptor heap for depth buffer.
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_DSVHeap;

    // Root signatures
    Microsoft::WRL::ComPtr<ID3D12RootSignature> m_RootSignature;

    // Pipeline state objects.
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_PipelineState;

    D3D12_VIEWPORT m_Viewport;
    D3D12_RECT m_ScissorRect;

    float m_FoV;

    DirectX::XMMATRIX m_ViewMatrix;
    DirectX::XMMATRIX m_ProjectionMatrix;

    Microsoft::WRL::ComPtr<ID3DBlob> m_pixelShaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> m_meshShaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> m_taskShaderBlob;


    bool        m_ContentLoaded         = false;

    float       m_ClearColor[4]         = {0.09f, 0.09f, 0.12f, 1.0f};
    double      m_fps                   = 0;
    double      m_frameTime             = 0;
    double      m_totalRunTime          = 0.0;
    bool        m_wireframe             = true;
    bool        m_frustumCulling        = true;
    ShadingMode m_shadingMode           = LOD_SHADING;
    float       m_LoDScale              = 1;
    float       m_debugFloatSlider      = 1;
    bool        m_geo_morphing          = true;
    bool        m_tessellation          = true;
    bool        m_screen_space_LoD      = false;
    bool        m_backFaceCulling       = true;
    bool        m_triplanarMapping      = false;
    bool        m_allowLighting         = true;
    float       m_triplanarScale        = 1;
    float       m_triplanarBlendGrade   = 5;
    float       m_displacementScale     = 0.0;

    // camera related variables
    float3      m_cameraPos             = float3(0, 0, 0);
    float       m_CameraYaw             = 0;
    float       m_CameraRoll            = 0;
    float       m_CameraSensibility     = 0.003f;
    float       m_cameraSpeed           = 20;
    bool        m_freeCamera            = false;
    bool        m_autoRotateScene       = true;
    float       m_autoCameraDistance    = 4.5f;
    float       m_CameraScrollScale     = 0.25f;
    float       m_autoRotationOffset    = 4;
    bool        m_lockCameraShaderConstant = false;
    float3      m_lockedCameraPos       = float3(0, 0, 0);
    
    // Descriptor heap
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_CBV_SRV_UAV_Heap;

    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_Sampler_Heap;


    // buffers
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_MeshletBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_VertexIndicesBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_PrimitiveIndicesBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_MorphIndicesBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_VertexBuffers;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_ObjectsBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_ConstantsBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_WorkQueueBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_WorkQueueCountersBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_WorkQueueCountersClearValuesBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource>              m_GlobalMeshPayloadBuffer;
    UINT8* m_mappedConstantData                         = nullptr;

    Microsoft::WRL::ComPtr<ID3D12Resource>              m_HeightMapTexture;


    // gpu handles
    D3D12_GPU_DESCRIPTOR_HANDLE m_MeshletsSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_VertexIndicesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_PrimitiveIndicesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_MorphIndicesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_VerticesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_ObjectsSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_WorkQueueSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_WorkQueueCountersSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_WorkQueueCountersClearValuesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_GlobalMeshPayloadSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_HeightMapTextureSrvHandle;

    Microsoft::WRL::ComPtr<ID3D12CommandSignature> m_commandSignature;

    Scene m_scene;
    char  m_model_file_path[1024] = "./assets/scenes/TorusInstances.glb";//"C:/Users/wolfr/Desktop/plane.glb";//"C:/Users/wolfr/Desktop/dragon.glb";//"./assets/scenes/TestScene.glb";
};