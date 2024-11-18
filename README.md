# meshlet-lod
## Setup Instructions
Install the following packages using vcpkg:<br>
```
vcpkg install imgui
vcpkg install directx-headers
vcpkg install meshoptimizer
vcpkg install assimp
```
Generate visual studio project files using cmake:<br>
```
cmake . -B build -A x64
```
