cmake_minimum_required( VERSION 3.10.1 ) # Latest version of CMake when this file was created.

# Set the C++ standard to C++17 (example)
set(CMAKE_CXX_STANDARD 17)

# Ensure the standard is enforced
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# copy assets to output directory
set(ASSET_DIR "${CMAKE_SOURCE_DIR}/MeshletLoD/assets")
set(OUTPUT_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>")


add_custom_command(
    OUTPUT ${OUTPUT_DIR}/assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${ASSET_DIR} ${OUTPUT_DIR}/assets
    COMMENT "Copying assets to output directory"
)

add_custom_target(CopyAssets ALL
    DEPENDS ${OUTPUT_DIR}/assets
)

# header files
set( HEADER_FILES
    inc/MeshletLoD.h
    inc/MeshletScene.h
    inc/MeshletMesh.h
    inc/HLSLnames.h
    shaders/Structures.fxh
    inc/Bone.h
    inc/Animation.h
    inc/Animator.h
    inc/Mesh.h
    inc/Scene.h
    inc/ViewDependentMeshletLoD.h
    shaders/ViewDependentStructures.fxh

)

# source files
set( SRC_FILES
    src/main.cpp
    src/MeshletLoD.cpp
    src/MeshletScene.cpp
    src/MeshletMesh.cpp
    src/Bone.cpp
    src/Animation.cpp
    src/Animator.cpp
    src/Mesh.cpp
    src/Scene.cpp
    src/ViewDependentMeshletLoD.cpp
)

# shaders
set( SHADER_FILES
    shaders/TaskShader.hlsl
    shaders/MeshShader.hlsl
    shaders/PixelShader.hlsl
    shaders/ObjectCulling_ComputeShader.hlsl
    shaders/ViewDependentMeshShader.hlsl
    shaders/ViewDependentTaskShader.hlsl
    shaders/ViewDependentPixelShader.hlsl
)

source_group( "Resources\\Shaders" FILES ${SHADER_FILES} )

set_source_files_properties( shaders/MeshShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Mesh
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/TaskShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Amplification
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/PixelShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Pixel
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/ObjectCulling_ComputeShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Compute
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/ViewDependentMeshShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Mesh
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/ViewDependentTaskShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Amplification
    VS_SHADER_MODEL 6.6
)

set_source_files_properties( shaders/ViewDependentPixelShader.hlsl PROPERTIES 
    VS_SHADER_TYPE Pixel
    VS_SHADER_MODEL 6.6
)

add_executable( MeshletLoD WIN32
    ${HEADER_FILES} 
    ${SRC_FILES}
    ${SHADER_FILES}
)

target_include_directories( MeshletLoD
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc
)

target_link_libraries( MeshletLoD
    DX12Lib
    D3DCompiler.lib
    Shlwapi.lib
)

add_dependencies(MeshletLoD CopyAssets)
