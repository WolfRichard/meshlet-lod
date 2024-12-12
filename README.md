# meshlet-lod
## Requirements
- DirectX 12 Ultimate compliant GPU
- NVIDIA RTX 20 series and above or AMD RDNA 2 architecture and newer
- Windows 10 / 11
- Visual Studio 2019 / 2022
- CMake & vcpgk
## Setup Instructions
### Install the following packages via vcpkg:<be>
```
vcpkg install imgui[core,dx12-binding,win32-binding]:x64-windows
vcpkg install directx-headers
vcpkg install meshoptimizer
vcpkg install assimp
```
### Generate visual studio project files via cmake:<be>
```
cmake . -B build -A x64
```
## Screen Shot
![failed to load preview](MeshletLoD/assets/screen_shot.png)
