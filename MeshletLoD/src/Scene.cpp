#include "Scene.h"
#include "Windows.h"
#include <random>
#if defined(max)
#undef max
#endif

using namespace DirectX;

Scene::Scene()
{
}

void Scene::init(std::string file_path)
{
    m_scene_objects.clear();
    for (auto mesh : m_meshes) delete mesh;
    m_meshes.clear();
    m_draw_task_count = 0;
    m_vertex_count = 0;
    m_triangles_count = 0;
    m_meshlet_counts.clear();
    m_indirect_attributes.clear();

    loadScene(file_path);
    
    
    
}


void Scene::processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform)
{
   
}


void Scene::loadScene(std::string file_path)
{
    
}

