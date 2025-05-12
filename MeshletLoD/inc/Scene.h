#pragma once

#include "Mesh.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include "../shaders/Structures.fxh"

#include <chrono>


// Loads and manages all unique meshes and their scene object instances
class Scene
{
public:
    Scene();
    
    // clears all allocated memory specific to the current scene, thus freeing up system memory as well as allowing to load a new scene
    void free();

    // frees and loads the specified scene file
    void init(std::string file_path);

    // scene object instances
    std::vector<S_SceneObject> m_scene_objects;
    
    // mesh data
    std::vector<std::vector<S_Vertex>>        m_vertices;
    std::vector<std::vector<uint>>            m_vertex_indices;
    std::vector<std::vector<unsigned char>>   m_primitive_indices;
    std::vector<std::vector<uint>>            m_morph_indices;
    std::vector<std::vector<S_Meshlet>>       m_meshlets;
 
    // scene meta data
    uint m_total_meshlet_count;
    uint m_totoal_triangle_count;
    uint m_total_vertex_count;
    uint m_mesh_count;

    // pre-compute benchmark statistics
    std::chrono::duration<double>              m_preProcessingTime;

private:
    // loads the specified scene file, will always look for a previously created back-up for scenes that have already have been loaded once
    // if no back-up file was found, it will pre-process the scene as usual and setup a new back-up file for the given scene
    void loadScene(std::string file_path);

    // helper function to recursivly itterate through the entire scene tree, to load all scene objects
    void processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform);

    // load and store a back-upfile to/from the current scene class
    void storeSceneToBackUp(std::string file_path);
    bool loadSceneFromBackUp(std::string file_path);

    // uinique meshes of the scene
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