# meshlet-lod
### setup instructions
install the following packages via vcpkg:
vcpkg install imgui
vcpkg install directx-headers
vcpkg install meshoptimizer
vcpkg install assimp

generate visual studio project files via cmake:
cmake . -B build -A x64
