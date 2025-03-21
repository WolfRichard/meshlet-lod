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
    std::vector<std::vector<unsigned int>*>    m_vertex_indices;
    std::vector<std::vector<unsigned char>*>   m_primitive_indices;
    std::vector<std::vector<S_Meshlet>*>       m_meshlets;
    std::vector<std::vector<S_MeshletGroup>*>  m_meshlet_groups;

private:
    void loadScene(std::string file_path);
    void processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform);
    std::vector<Mesh*> m_meshes;
};

