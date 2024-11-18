# meshlet-lod
## setup instructions
###install the following packages via vcpkg:<br>
vcpkg install imgui<br>
vcpkg install directx-headers<br>
vcpkg install meshoptimizer<br>
vcpkg install assimp<br>

###generate visual studio project files via cmake:<br>
cmake . -B build -A x64<br>
