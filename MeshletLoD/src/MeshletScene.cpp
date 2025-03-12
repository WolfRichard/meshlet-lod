#include "MeshletScene.h"
#include "Windows.h"
#include <random>
#if defined(max)
#undef max
#endif

using namespace DirectX;

MeshletScene::MeshletScene()
{
}

void MeshletScene::init(std::string file_path, uint selectedLoD)
{
    m_scene_objects.clear();
    for (auto mesh : m_meshes) delete mesh;
    m_meshes.clear();
    m_draw_task_count = 0;
    m_vertex_count = 0;
    m_triangles_count = 0;
    m_mesh_lod_count = 0;
    m_meshlet_counts.clear();
    m_indirect_attributes.clear();
    m_preBakedAnimations.clear();
    m_animationMetaData.clear();
    m_sceneObjectNames.clear();
    m_sceneObjectNamesCharP.clear();
    m_meshLoDBufferStructure.clear();

    m_modelLoadTime = m_modelLoadTime.zero();
    m_totalLoDGenTime = m_totalLoDGenTime.zero();
    m_totalMeshletGenTime = m_totalMeshletGenTime.zero();
    m_totalTime = m_totalTime.zero();

    auto benchmark_time_start = std::chrono::high_resolution_clock::now();
    loadScene(file_path, selectedLoD);
    
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_totalTime = benchmark_time_end - benchmark_time_start;
    
}


void MeshletScene::processSceneNode(aiNode* node, const aiScene* scene, float4x4 parent_transform)
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

        // animation data
        so.animation_id = -1;
        so.animation_speed = 1;
        so.animation_time_offset = 0;
        if (m_meshes[node->mMeshes[i]]->m_animations.size() > 0)
            so.animation_id = m_meshes[node->mMeshes[i]]->m_animations[0].totalAnimationIndex;

        m_scene_objects.push_back(so);
        m_draw_task_count += m_meshlet_counts[m_meshLoDBufferStructure[node->mMeshes[i]].mesh_offset];
        m_vertex_count += (uint)m_meshes[node->mMeshes[i]]->m_vertex_count;
        m_triangles_count += (uint)m_meshes[node->mMeshes[i]]->m_indices.size() / 3;

        uint thread_count = ((m_meshlet_counts[node->mMeshes[i]] + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE;

        m_indirect_attributes.push_back(CommandStructure());
        m_indirect_attributes.back().instanceID = (uint)(m_scene_objects.size() - 1u);
        m_indirect_attributes.back().level_of_detail = 0;
        m_indirect_attributes.back().dispatchArguments.x = thread_count;
        m_indirect_attributes.back().dispatchArguments.y = 1u;
        m_indirect_attributes.back().dispatchArguments.z = 1u;

        
        std::string nodeName = std::string("Object_" + std::to_string(m_scene_objects.size() - 1));
        if (node->mName.length)
            nodeName = std::string(node->mName.C_Str());
        m_sceneObjectNames.push_back(nodeName);
        m_sceneObjectNamesCharP.push_back(m_sceneObjectNames.back().c_str());


        //genInstances((int)m_scene_objects.size() - 1, 5, 5);
    }

    for (uint i = 0; i < node->mNumChildren; i++)
    {
        processSceneNode(node->mChildren[i], scene, model_matrix);
    }
}


void MeshletScene::loadScene(std::string file_path, uint selectedLoD)
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
            m_meshes.back()->m_animations[a].totalAnimationIndex = (uint)m_preBakedAnimations.size();
            m_preBakedAnimations.push_back(&(m_meshes.back()->m_animations[a]));
            AnimationMetaData amd;
            amd.bone_count = m_meshes.back()->m_animations[a].boneCount;
            amd.frame_count = m_meshes.back()->m_animations[a].frameCount;
            amd.duration = (float)(m_meshes.back()->m_animations[a].frameCount - 1) / (float)ANIMATION_FPS;
            m_animationMetaData.push_back(amd);
        }

        m_mesh_lod_count += m_meshes.back()->m_LoDCount;

        MeshLoDStructure mls;
        mls.lod_count = m_meshes.back()->m_LoDCount;
        mls.mesh_offset = (uint)m_meshlet_counts.size();
        m_meshLoDBufferStructure.push_back(mls);

        for (uint i = 0; i < (uint)m_meshes.back()->m_LoDCount; i++)
        {
            m_meshlet_counts.push_back((uint)m_meshes.back()->m_draw_tasks[i].size());
        }

        m_totalLoDGenTime += m_meshes.back()->m_LoDGenTime;
        m_totalMeshletGenTime += m_meshes.back()->m_meshletGenTime;
    }


    // process scene tree
    float4x4 base_transform = float4x4(1, 0, 0, 0,
                                       0, 1, 0, 0,
                                       0, 0, 1, 0,
                                       0, 0, 0, 1);
    processSceneNode(scene->mRootNode, scene, base_transform);
}

float GetRandomFloat(float min, float max) {
    return min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min)));
}

int GetRandomInt(int min, int max) {
    static std::random_device rd;
    static std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

void MeshletScene::genInstances(uint objectIndex, uint gridSideCount, float spacing)
{
    SceneObject selectedObject = m_scene_objects[objectIndex];
    std::string objectName = m_sceneObjectNames[objectIndex];

    for (uint x = 0; x < gridSideCount; x++)
    {
        for (uint z = 0; z < gridSideCount; z++)
        {
            if (!x && !z)
                continue;
            SceneObject so;
            so.object_matrix = DirectX::XMMatrixTranspose(DirectX::XMMatrixTranslation(x * spacing, 0, z * spacing)) * selectedObject.object_matrix;
            so.mesh_id = selectedObject.mesh_id;
            so.bounding_sphere_center = selectedObject.bounding_sphere_center;
            so.bounding_sphere_radius = selectedObject.bounding_sphere_radius;

            // animation data
            so.animation_id = selectedObject.animation_id;
            so.animation_speed = selectedObject.animation_speed;
            so.animation_time_offset = selectedObject.animation_time_offset;

            m_scene_objects.push_back(so);
            m_draw_task_count += m_meshlet_counts[m_meshLoDBufferStructure[so.mesh_id].mesh_offset];
            m_vertex_count += (uint)m_meshes[so.mesh_id]->m_vertex_count;
            m_triangles_count += (uint)m_meshes[so.mesh_id]->m_meshlet_triangles.size() / 3;

            uint thread_count = ((m_meshlet_counts[so.mesh_id] + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE;

            m_indirect_attributes.push_back(CommandStructure());
            m_indirect_attributes.back().instanceID = (uint)(m_scene_objects.size() - 1u);
            m_indirect_attributes.back().level_of_detail = 0;
            m_indirect_attributes.back().dispatchArguments.x = thread_count;
            m_indirect_attributes.back().dispatchArguments.y = 1u;
            m_indirect_attributes.back().dispatchArguments.z = 1u;

            m_sceneObjectNames.push_back(std::string("Instance_" + std::to_string(x * gridSideCount + z)));
            m_sceneObjectNamesCharP.push_back(m_sceneObjectNames.back().c_str());
        }
    }
}

void MeshletScene::randomizeAnimationOffsets()
{
    for (uint o  = 0; o < (uint)m_scene_objects.size(); o++)
        m_scene_objects[o].animation_time_offset = GetRandomFloat(0, 10.0f);
}

void MeshletScene::randomizeAnimations()
{
    for (uint o = 0; o < (uint)m_scene_objects.size(); o++)
    {
        if (m_meshes[m_scene_objects[o].mesh_id]->m_animations.size())
        {
            m_scene_objects[o].animation_id = m_meshes[m_scene_objects[o].mesh_id]->m_animations[GetRandomInt(0, (int)m_meshes[m_scene_objects[o].mesh_id]->m_animations.size() - 1)].totalAnimationIndex;
            m_scene_objects[o].animation_time_offset = GetRandomFloat(0, 10.0f);
            m_scene_objects[o].animation_speed = GetRandomFloat(0.5f, 2.0f);
        }
    }

}

void MeshletScene::resetAnimations()
{
    for (uint o = 0; o < (uint)m_scene_objects.size(); o++)
    {
        if (m_meshes[m_scene_objects[o].mesh_id]->m_animations.size())
        {
            m_scene_objects[o].animation_id = m_meshes[m_scene_objects[o].mesh_id]->m_animations[0].totalAnimationIndex;
            m_scene_objects[o].animation_time_offset = 0;
            m_scene_objects[o].animation_speed = 1;
        }
    }
}

void MeshletScene::clearAnimations()
{
    for (uint o = 0; o < (uint)m_scene_objects.size(); o++)
        m_scene_objects[o].animation_id = -1;
}

