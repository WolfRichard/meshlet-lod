#include "Scene.h"

using namespace DirectX;

Scene::Scene()
{
}

void Scene::init(std::string file_path, uint selectedLoD)
{
    loadScene(file_path, selectedLoD);
}


void Scene::processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform)
{
    aiMatrix4x4 t              = node->mTransformation;
    float4x4    node_transform = float4x4(t.a1, t.a2, t.a3, t.a4,
                                          t.b1, t.b2, t.b3, t.b4,
                                          t.c1, t.c2, t.c3, t.c4,
                                          t.d1, t.d2, t.d3, t.d4);

    float4x4 model_matrix = node_transform * parent_transform;

    for (uint i = 0; i < node->mNumMeshes; i++)
    {
        SceneObject so;
        so.object_matrix = model_matrix;
        so.mesh_id       = node->mMeshes[i];

        // transform object space bounding sphere ti world space

        XMVECTOR wSpaceCentre = XMVector3Transform(XMLoadFloat3(&(m_meshes[so.mesh_id]->m_boundingSphereCentre)), model_matrix);
        XMStoreFloat3(&so.bounding_sphere_center, wSpaceCentre);
        // scale bounding radius by highest axis scale factor
        XMVECTOR scaleX = XMVector3Length(model_matrix.r[0]); 
        XMVECTOR scaleY = XMVector3Length(model_matrix.r[1]); 
        XMVECTOR scaleZ = XMVector3Length(model_matrix.r[2]); 
        float maxScale = std::max({ XMVectorGetX(scaleX), XMVectorGetX(scaleY), XMVectorGetX(scaleZ) });
        so.bounding_sphere_radius = m_meshes[so.mesh_id]->m_boundingSphereRadius * maxScale;

        m_scene_objects.push_back(so);
        m_draw_task_count += m_meshlet_counts[node->mMeshes[i]];

        uint thread_count = ((m_meshlet_counts[node->mMeshes[i]] + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE;
        m_indirect_attributes.push_back(thread_count);
        m_indirect_attributes.push_back(1);
        m_indirect_attributes.push_back(1);

        m_indirect_attributes_withConstant.push_back((uint)m_scene_objects.size() - 1);
        m_indirect_attributes_withConstant.push_back(thread_count);
        m_indirect_attributes_withConstant.push_back(1);
        m_indirect_attributes_withConstant.push_back(1);
        
    }

    for (uint i = 0; i < node->mNumChildren; i++)
    {
        processSceneNode(node->mChildren[i], scene, model_matrix);
    }
}


void Scene::loadScene(std::string file_path, uint selectedLoD)
{
    // try to load the model via assimp
    Assimp::Importer importer;
    const aiScene*   scene = importer.ReadFile(file_path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenSmoothNormals); // | aiProcess_RemoveRedundantMaterials | aiProcess_OptimizeMeshes | aiProcess_OptimizeGraph | aiProcess_JoinIdenticalVertices);

    // verify that loading was successfull + error handeling
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        throw std::exception();
    }


    // load meshes
    for (uint m = 0; m < scene->mNumMeshes; m++)
    {
        m_meshes.push_back(new MeshletMesh(scene->mMeshes[m]));

        m_meshlet_counts.push_back((uint)m_meshes.back()->m_draw_tasks.size());
    }

    // process scene tree
    float4x4 base_transform = float4x4(1, 0, 0, 0,
                                       0, 1, 0, 0,
                                       0, 0, 1, 0,
                                       0, 0, 0, 1);
    processSceneNode(scene->mRootNode, scene, base_transform);
}

