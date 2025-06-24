#pragma once

#include <Game.h>
#include <Window.h>
#include <vector>

#include <DirectXMath.h>

#include "MeshletScene.h"




class DiscreteMeshletLoD : public Game
{
public:
    using super = Game;

    DiscreteMeshletLoD(const std::wstring& name, int width, int height, bool vSync = false);
    /**
     *  Load content required for the demo.
     */
    virtual bool LoadContent() override;

    /**
     *  Unload demo specific content that was loaded in LoadContent.
     */
    virtual void UnloadContent() override;
protected:
    /**
     *  Update the game logic.
     */
    virtual void OnUpdate(UpdateEventArgs& e) override;

    /**
     *  Render stuff.
     */
    virtual void OnRender(RenderEventArgs& e) override;

    /**
     * Invoked by the registered window when a key is pressed
     * while the window has focus.
     */
    virtual void OnKeyPressed(KeyEventArgs& e) override;

    /**
     * Invoked when the mouse wheel is scrolled while the registered window has focus.
     */
    virtual void OnMouseWheel(MouseWheelEventArgs& e) override;


    virtual void OnResize(ResizeEventArgs& e) override; 

private:
    enum DebugVisualsSelection {
        ShowNoDebug,
        ShowMeshlets,
        ShowBones,
        ShowLoD
    };

    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_ImGuiDescriptorHeap;
    void initImGui();
    void updateImGui();
    void CreatePSO();
    void CreateCullingPSO();
    void UpdateObjectsBuffer();

    // Helper functions
    // Transition a resource
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
    Microsoft::WRL::ComPtr<ID3D12RootSignature> m_ObjectCulling_RootSignature;

    // Pipeline state objects.
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_PipelineState;
    Microsoft::WRL::ComPtr<ID3D12PipelineState> m_ObjectCulling_PipelineState;


    D3D12_VIEWPORT m_Viewport;
    D3D12_RECT m_ScissorRect;

    float m_FoV;

    DirectX::XMMATRIX m_ViewMatrix;
    DirectX::XMMATRIX m_ProjectionMatrix;

    Microsoft::WRL::ComPtr<ID3DBlob> m_pixelShaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> m_meshShaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> m_taskShaderBlob;
    Microsoft::WRL::ComPtr<ID3DBlob> m_objectCullingComputeShaderBlob;

    bool m_ContentLoaded;

    float m_ClearColor[4] = {23.0f / 255.0f, 23.0f / 255.0f, 31.0f / 255.0f, 1.0f};
    double m_fps = 0;
    double m_frameTime = 0;
    double m_totalRunTime = 0.0;
    bool m_wireframe = false;
    bool m_backFaceCulling = true;
    bool m_frustumCulling = true;
    bool m_coneCulling = false;
    DebugVisualsSelection m_debugMode = ShowNoDebug;
    bool m_objectCulling = true;
    bool m_objectLoD = true;
    bool m_meshletLoD = true;

    // camera related variables
    float3 m_cameraPos = float3(0, 0, 0);
    float       m_CameraYaw = 0;
    float       m_CameraRoll = 0;

    float m_CameraSensibility = 0.003f;
    float m_cameraSpeed = 20;

    
    bool m_freeCamera = false;
    bool m_autoRotateScene = true;
    float m_autoCameraDistance = 4.5f;
    float m_CameraScrollScale = 0.25f;
    float m_autoRotationOffset = 4;
    float m_LoDScale = 2;
    
    

    


    // Descriptor heap
    Microsoft::WRL::ComPtr<ID3D12DescriptorHeap> m_CBV_SRV_UAV_Heap;

    // buffers
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_IndexBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_VertexBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_TriangleBuffers;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_DrawTasksBuffers;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_ObjectsBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_MeshletCountsBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_indirectArgumentBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_objectCountBuffer;
    std::vector<Microsoft::WRL::ComPtr<ID3D12Resource>> m_BoneMatricesBuffers;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_AnimationMetaDataBuffer;
    Microsoft::WRL::ComPtr<ID3D12Resource> m_MeshLoDStructureBuffer;
    
    
    Microsoft::WRL::ComPtr<ID3D12Resource> m_ConstantsBuffer;
    UINT8* m_mappedConstantData = nullptr;


    D3D12_CPU_DESCRIPTOR_HANDLE objectBufferSRVHandle;

    // gpu handles
    D3D12_GPU_DESCRIPTOR_HANDLE m_indexSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_vertexSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_triangleSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_drawTasksSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_objectsSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_meshletCountsSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_visibleObjectCountSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_boneMatricesSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_animationMetaDataSrvHandle;
    D3D12_GPU_DESCRIPTOR_HANDLE m_meshLoDStructureSrvHandle;


    Microsoft::WRL::ComPtr<ID3D12CommandSignature> m_commandSignature;

    MeshletScene m_scene;
    char  m_model_file_path[1024] = "./assets/scenes/TestScene.glb";//"./assets/scenes/TorusInstances.glb";
};

//"C:/Users/wolfr/Downloads/xyzrgb_dragon.ply/xyzrgb_dragon.ply";