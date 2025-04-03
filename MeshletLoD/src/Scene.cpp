#include "Scene.h"
#include "Windows.h"
#include <random>
#if defined(max)
#undef max
#endif

using namespace DirectX;

Scene::Scene()
{
    m_mesh_count            = 0;
    m_total_meshlet_count   = 0;
    m_totoal_triangle_count = 0;
    m_total_vertex_count    = 0;
    m_preProcessingTime     = m_preProcessingTime.zero();
}

void Scene::init(std::string file_path)
{

    

    
    m_scene_objects.clear();
    for (auto mesh : m_meshes) delete mesh;
    m_meshes.clear();

    auto benchmark_time_start = std::chrono::high_resolution_clock::now();
    loadScene(file_path);
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_preProcessingTime = benchmark_time_end - benchmark_time_start;
}


void Scene::processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform)
{
    aiMatrix4x4 t = node->mTransformation;
    float4x4    node_transform = float4x4(t.a1, t.a2, t.a3, t.a4,
        t.b1, t.b2, t.b3, t.b4,
        t.c1, t.c2, t.c3, t.c4,
        t.d1, t.d2, t.d3, t.d4);

    float4x4 model_matrix = node_transform * parent_transform;

    for (uint i = 0; i < node->mNumMeshes; i++)
    {
        S_SceneObject so;
        so.object_matrix = model_matrix;
        so.mesh_id = node->mMeshes[i];

        // transform object space bounding sphere ti world space

        XMVECTOR wSpaceCentre = XMVector3Transform(XMLoadFloat3(&(m_meshes[so.mesh_id]->m_bounding_sphere.center)), XMMatrixTranspose(model_matrix));
        XMStoreFloat3(&so.bounding_sphere.center, wSpaceCentre);
        // scale bounding radius by highest axis scale factor
        XMVECTOR scaleX = XMVector3Length(model_matrix.r[0]);
        XMVECTOR scaleY = XMVector3Length(model_matrix.r[1]);
        XMVECTOR scaleZ = XMVector3Length(model_matrix.r[2]);
        float maxScale = std::max({ XMVectorGetX(scaleX), XMVectorGetX(scaleY), XMVectorGetX(scaleZ) });
        so.bounding_sphere.radius = m_meshes[so.mesh_id]->m_bounding_sphere.radius * maxScale;
        so.root_group_id = m_meshes[so.mesh_id]->m_hierarchy_root_group;
        so.mesh_group_count = (uint)m_meshes[so.mesh_id]->m_meshlet_groups.size();
        so.mesh_meshlet_count = (uint)m_meshes[so.mesh_id]->m_meshlets.size();

        
        m_total_meshlet_count   += (uint)m_meshes[so.mesh_id]->m_meshlets.size();
        m_totoal_triangle_count += (uint)m_meshes[so.mesh_id]->m_original_indices.size() / 3;
        m_total_vertex_count    += (uint)m_meshes[so.mesh_id]->m_vertices.size();



        m_scene_objects.push_back(so);
    }

    for (uint i = 0; i < node->mNumChildren; i++)
    {
        processSceneNode(node->mChildren[i], scene, model_matrix);
    }
}


void Scene::loadScene(std::string file_path)
{
    // try to load the model via assimp
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(file_path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenSmoothNormals | aiProcess_MakeLeftHanded | aiProcess_FlipWindingOrder | aiProcess_RemoveRedundantMaterials | aiProcess_OptimizeMeshes | aiProcess_OptimizeGraph | aiProcess_JoinIdenticalVertices);

    // verify that loading was successfull + error handeling
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        throw std::exception();
    }

    // load meshes
    for (uint m = 0; m < scene->mNumMeshes; m++)
    {
        m_meshes.push_back(new Mesh(scene->mMeshes[m], scene));

        m_vertices.push_back(&(m_meshes.back()->m_vertices));
        m_vertex_indices.push_back(&(m_meshes.back()->m_vertex_indices));
        m_primitive_indices.push_back(&(m_meshes.back()->m_primitive_indices));
        m_meshlets.push_back(&(m_meshes.back()->m_meshlets));
        m_meshlet_groups.push_back(&(m_meshes.back()->m_meshlet_groups));
        m_morph_indices.push_back(&(m_meshes.back())->m_morph_indices);
    }
    m_mesh_count = scene->mNumMeshes;

    // process scene tree
    float4x4 base_transform = float4x4(1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1);
    processSceneNode(scene->mRootNode, scene, base_transform);
}

