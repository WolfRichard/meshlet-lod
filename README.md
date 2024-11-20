# meshlet-lod
![failed to load preview](MeshletLoD/assets/screen_shot.png)

## requirements
- DirectX 12 Ultimate compliant GPU
- NVIDIA RTX 20 series and above or AMD RDNA 2 architecture and newer
- Windows 10 / 11
- Visual Studio 2019 / 2022

## setup instructions
### install the following packages via vcpkg:<br>
vcpkg install imgui<br>
vcpkg install directx-headers<br>
vcpkg install meshoptimizer<br>
vcpkg install assimp<br>

### generate visual studio project files via cmake:<br>
cmake . -B build -A x64<br>
