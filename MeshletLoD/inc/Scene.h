#pragma once

#include "Mesh.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include "../shaders/ViewDependentStructures.fxh"

#include <chrono>


class Scene
{
public:
    Scene();

    void init(std::string file_path);
    
    uint m_mesh_count;

    std::vector<S_SceneObject> m_scene_objects;
    
    std::vector<std::vector<S_Vertex>*>        m_vertices;
    std::vector<std::vector<uint>*>            m_vertex_indices;
    std::vector<std::vector<unsigned char>*>   m_primitive_indices;
    std::vector<std::vector<uint>*>            m_morph_indices;
    std::vector<std::vector<S_Meshlet>*>       m_meshlets;
    //std::vector<std::vector<S_MeshletGroup>*>  m_meshlet_groups;


    // statistics
    uint m_total_meshlet_count;
    uint m_totoal_triangle_count;
    uint m_total_vertex_count;
    std::chrono::duration<double>              m_preProcessingTime;

private:
    void loadScene(std::string file_path);
    void processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform);

    void storeSceneToBackUp(std::string file_path);
    bool loadSceneFromBackUp(std::string file_path);

    std::vector<Mesh*> m_meshes;
};

/******************************************************************
                Pre-Processing Back-Up File Layout
*******************************************************************

uint32 scene_object_count;
uint32 unique_mesh_count;

S_SceneObject scene_objects_buffer[scene_object_count];

for_each (unique mesh) {
    uint vertex_count
    S_Vertex m_vertices[vertex_count]
}

for_each (unique mesh) {
    uint vertex_indices_count
    uint m_vertex_indices[vertex_indices_count]
}

for_each (unique mesh) {
    uint morph_indices_count
    uint m_morph_indices[morph_indices_count]
}

for_each (unique mesh) {
    uint primitive_indices_count
    unsigned char m_primitive_indices[primitive_indices_count]
}

for_each (unique mesh) {
    uint meshlet_count
    S_Meshlet m_meshlets[meshlet_count]
}
*/