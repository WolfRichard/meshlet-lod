#pragma once

#include "MeshletMesh.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include "../shaders/Structures.fxh"



class Scene
{
public:
    Scene();

    void init(std::string file_path, uint selectedLoD);

    std::vector<SceneObject> m_scene_objects;
    std::vector<MeshletMesh*> m_meshes;

    uint m_draw_task_count = 0;
    uint m_vertex_count = 0;
    uint m_triangles_count = 0;

    std::vector<uint> m_meshlet_counts;
    std::vector<CommandStructure> m_indirect_attributes;

private:
    void loadScene(std::string file_path, uint selectedLoD);
    void processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform);
};

