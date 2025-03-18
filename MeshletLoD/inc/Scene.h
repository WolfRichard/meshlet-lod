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
    
    std::vector<S_SceneObject> m_scene_objects;
    std::vector<Mesh*> m_meshes;
   

private:
    void loadScene(std::string file_path);
    void processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform);
};

