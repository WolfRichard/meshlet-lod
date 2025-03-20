#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <Shlwapi.h>

#include <Application.h>
#include <MeshletLoD.h>
#include <ViewDependentMeshletLoD.h>

#include <dxgidebug.h>

void ReportLiveObjects()
{
    IDXGIDebug1* dxgiDebug;
    DXGIGetDebugInterface1(0, IID_PPV_ARGS(&dxgiDebug));

    dxgiDebug->ReportLiveObjects(DXGI_DEBUG_ALL, DXGI_DEBUG_RLO_IGNORE_INTERNAL);
    dxgiDebug->Release();
}

int CALLBACK wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR lpCmdLine, int nCmdShow)
{
    int retCode = 0;

    // Set the working directory to the path of the executable.
    WCHAR path[MAX_PATH];
    HMODULE hModule = GetModuleHandleW(NULL);
    if ( GetModuleFileNameW(hModule, path, MAX_PATH) > 0 )
    {
        PathRemoveFileSpecW(path);
        SetCurrentDirectoryW(path);
    }

    Application::Create(hInstance);
    {
        std::shared_ptr<MeshletLoD> game = std::make_shared<MeshletLoD>(L"Meshlet_LoD", 1280, 720);
        //std::shared_ptr<ViewDependentMeshletLoD> game = std::make_shared<ViewDependentMeshletLoD>(L"ViewDependentMeshlet_LoD", 1280, 720);
        retCode = Application::Get().Run(game);
    }
    Application::Destroy();

    atexit(&ReportLiveObjects);

    return retCode;
}