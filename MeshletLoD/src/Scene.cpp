#include "Scene.h"

using namespace DirectX;

Scene::Scene()
{
}

void Scene::init(std::string file_path, uint selectedLoD)
{
    m_scene_objects.clear();
    for (auto mesh : m_meshes) delete mesh;
    m_meshes.clear();
    m_draw_task_count = 0;
    m_vertex_count = 0;
    m_triangles_count = 0;
    m_meshlet_counts.clear();
    m_indirect_attributes.clear();

    m_modelLoadTime = m_modelLoadTime.zero();
    m_totalLoDGenTime = m_totalLoDGenTime.zero();
    m_totalMeshletGenTime = m_totalMeshletGenTime.zero();
    m_totalTime = m_totalTime.zero();

    auto benchmark_time_start = std::chrono::high_resolution_clock::now();
    loadScene(file_path, selectedLoD);
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_totalTime = benchmark_time_end - benchmark_time_start;
    
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

        XMVECTOR wSpaceCentre = XMVector3Transform(XMLoadFloat3(&(m_meshes[so.mesh_id]->m_boundingSphereCentre)), XMMatrixTranspose(model_matrix));
        XMStoreFloat3(&so.bounding_sphere_center, wSpaceCentre);
        // scale bounding radius by highest axis scale factor
        XMVECTOR scaleX = XMVector3Length(model_matrix.r[0]); 
        XMVECTOR scaleY = XMVector3Length(model_matrix.r[1]); 
        XMVECTOR scaleZ = XMVector3Length(model_matrix.r[2]); 
        float maxScale = std::max({ XMVectorGetX(scaleX), XMVectorGetX(scaleY), XMVectorGetX(scaleZ) });
        so.bounding_sphere_radius = m_meshes[so.mesh_id]->m_boundingSphereRadius * maxScale;

        m_scene_objects.push_back(so);
        m_draw_task_count += m_meshlet_counts[node->mMeshes[i]];
        m_vertex_count += (uint)m_meshes[node->mMeshes[i]]->m_vertex_count;
        m_triangles_count += (uint)m_meshes[node->mMeshes[i]]->m_index_count / 3;

        uint thread_count = ((m_meshlet_counts[node->mMeshes[i]] + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE;

        m_indirect_attributes.push_back(CommandStructure());
        m_indirect_attributes.back().instanceID = (uint)(m_scene_objects.size() - 1u);
        m_indirect_attributes.back().dispatchArguments.x = thread_count;
        m_indirect_attributes.back().dispatchArguments.y = 1u;
        m_indirect_attributes.back().dispatchArguments.z = 1u;

        // animation data
        so.animation_id = -1;
        so.animation_speed = 1;
        so.animation_time_offset = 0;

        if (m_meshes[node->mMeshes[i]]->m_animations.size() > 0)
        {
            so.animation_id = m_meshes[node->mMeshes[i]]->m_animations[0].totalAnimationIndex;

        }
    }

    for (uint i = 0; i < node->mNumChildren; i++)
    {
        processSceneNode(node->mChildren[i], scene, model_matrix);
    }
}


void Scene::loadScene(std::string file_path, uint selectedLoD)
{
    // try to load the model via assimp
    auto benchmark_time_start = std::chrono::high_resolution_clock::now();
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(file_path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenSmoothNormals | aiProcess_MakeLeftHanded | aiProcess_FlipWindingOrder); // | aiProcess_RemoveRedundantMaterials | aiProcess_OptimizeMeshes | aiProcess_OptimizeGraph | aiProcess_JoinIdenticalVertices);

    // verify that loading was successfull + error handeling
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        throw std::exception();
    }
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();

    m_modelLoadTime = benchmark_time_end - benchmark_time_start;

    // load meshes
    for (uint m = 0; m < scene->mNumMeshes; m++)
    {
        m_meshes.push_back(new MeshletMesh(scene->mMeshes[m], scene));

        // animation 
        for (auto a = 0; a < m_meshes.back()->m_animations.size(); a++)
        {
            m_meshes.back()->m_animations[a].totalAnimationIndex = m_preBakedAnimations.size();
            m_preBakedAnimations.push_back(&(m_meshes.back()->m_animations[a]));
            AnimationMetaData amd;
            amd.bone_count = m_meshes.back()->m_animations[a].boneCount;
            amd.frame_count = m_meshes.back()->m_animations[a].frameCount;
            amd.duration = (float)(m_meshes.back()->m_animations[a].frameCount - 1) / (float)ANIMATION_FPS;
            m_animationMetaData.push_back(amd);
        }

        m_meshlet_counts.push_back((uint)m_meshes.back()->m_draw_tasks.size());

        m_totalLoDGenTime += m_meshes.back()->m_LoDGenTime;
        m_totalMeshletGenTime += m_meshes.back()->m_meshletGenTime;
    }

    if (scene->HasAnimations())
    {
        testAnimation.init(scene, 0, m_meshes[0]);
        animator.init(&testAnimation);
    }
    else
        animator.init(nullptr);
    //animator.PlayAnimation(&testAnimation);

    // process scene tree
    float4x4 base_transform = float4x4(1, 0, 0, 0,
                                       0, 1, 0, 0,
                                       0, 0, 1, 0,
                                       0, 0, 0, 1);
    processSceneNode(scene->mRootNode, scene, base_transform);
}


