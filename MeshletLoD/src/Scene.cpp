#include "Scene.h"
#include "Windows.h"
#include <filesystem>
#include <fstream>
#include <random>
#if defined(max)
#undef max
#endif

using namespace DirectX;
namespace fs = std::filesystem;

Scene::Scene()
{
    free();
}

void Scene::init(std::string file_path)
{
    free();

    auto benchmark_time_start = std::chrono::high_resolution_clock::now();


    // generate backup file name from original file path
    std::string extension = fs::path(file_path).extension().string();
    std::string filename_without_ext = file_path.substr(0, file_path.size() - extension.size());
    std::string backup_path = filename_without_ext + ".lod_bak";

    // check if backup file for precomuted scene can be found, if so load scene directly from backup and skip preprocessing
    if (!loadSceneFromBackUp(backup_path))
    {
        // preprocess and store data inside a backup file for future reuse
        loadScene(file_path);
        storeSceneToBackUp(backup_path);
    }

    
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_totalPreProcessingTime = benchmark_time_end - benchmark_time_start;
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
        so.mesh_meshlet_count = (uint)m_meshes[so.mesh_id]->m_meshlets.size();

        // update scene meta data counters
        m_total_meshlet_count   += (uint)m_meshes[so.mesh_id]->m_meshlets.size();
        m_totoal_triangle_count += (uint)m_meshes[so.mesh_id]->m_original_indices.size() / 3;
        m_total_vertex_count    += (uint)m_meshes[so.mesh_id]->m_vertices.size();
        
        // store the fully loaded scene object instance
        m_scene_objects.push_back(so);
    }

    // repeat the process for each child node of the currently processed node until whole tree is traversed
    for (uint i = 0; i < node->mNumChildren; i++)
    {
        processSceneNode(node->mChildren[i], scene, model_matrix);
    }
}


void Scene::loadScene(std::string file_path)
{
    auto scene_import_time_start = std::chrono::high_resolution_clock::now();

    // try to load the model via assimp
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(file_path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_ForceGenNormals | aiProcess_GenSmoothNormals | aiProcess_MakeLeftHanded | aiProcess_FlipWindingOrder | aiProcess_RemoveRedundantMaterials | aiProcess_JoinIdenticalVertices);

    // verify that loading was successfull + error handeling
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
    {
        throw std::exception();
    }

    auto scene_import_time_end = std::chrono::high_resolution_clock::now();
    m_sceneImportTime = scene_import_time_end - scene_import_time_start;

    m_modelParsingTime = m_modelParsingTime.zero();
    m_hierarchyGenTime = m_hierarchyGenTime.zero();
    
    // load meshes
    for (uint m = 0; m < scene->mNumMeshes; m++)
    {
        m_meshes.push_back(new Mesh(scene->mMeshes[m], scene));

        m_vertices.push_back(m_meshes.back()->m_vertices);
        m_vertex_indices.push_back(m_meshes.back()->m_vertex_indices);
        m_primitive_indices.push_back(m_meshes.back()->m_primitive_indices);
        m_meshlets.push_back(m_meshes.back()->m_meshlets);
        m_morph_indices.push_back(m_meshes.back()->m_morph_indices);

        m_modelParsingTime += m_meshes.back()->m_meshParseTime;
        m_hierarchyGenTime += m_meshes.back()->m_hierarchyGenTime;
    }
    m_mesh_count = scene->mNumMeshes;

 
    


    // process scene tree
    scene_import_time_start = std::chrono::high_resolution_clock::now();
    float4x4 base_transform = float4x4(1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1);
    processSceneNode(scene->mRootNode, scene, base_transform);

    scene_import_time_end = std::chrono::high_resolution_clock::now();
    m_sceneImportTime += scene_import_time_end - scene_import_time_start;
}


void Scene::storeSceneToBackUp(std::string file_path)
{
    std::ofstream back_up_file(file_path, std::ios::binary);
    if (!back_up_file)
    {
        OutputDebugString("FAILED TO OPEN/CREATE BACKUP FILE FOR LATER RE-USE!\n");
        return;
    }

    // backup scene meta data
    uint scene_object_count = (uint)m_scene_objects.size();
    back_up_file.write(reinterpret_cast<const char*>(&scene_object_count), sizeof(uint));
    back_up_file.write(reinterpret_cast<const char*>(&m_mesh_count), sizeof(uint));

    // backup scene objects buffer data
    back_up_file.write(reinterpret_cast<const char*>(m_scene_objects.data()), sizeof(S_SceneObject) * m_scene_objects.size());
    
    // backup vertex buffers
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint vertex_count = (uint)m_vertices[i].size();
        back_up_file.write(reinterpret_cast<const char*>(&vertex_count), sizeof(uint));
        back_up_file.write(reinterpret_cast<const char*>(m_vertices[i].data()), sizeof(S_Vertex) * vertex_count);
    }

    // backup vertex-indices buffers
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint vertex_indices_count = (uint)m_vertex_indices[i].size();
        back_up_file.write(reinterpret_cast<const char*>(&vertex_indices_count), sizeof(uint));
        back_up_file.write(reinterpret_cast<const char*>(m_vertex_indices[i].data()), sizeof(uint) * vertex_indices_count);
    }

    // backup morph-indices buffers
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint morph_indices_count = (uint)m_morph_indices[i].size();
        back_up_file.write(reinterpret_cast<const char*>(&morph_indices_count), sizeof(uint));
        back_up_file.write(reinterpret_cast<const char*>(m_morph_indices[i].data()), sizeof(uint) * morph_indices_count);
    }

    // backup primitive_indices buffers
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint primitive_indices_count = (uint)m_primitive_indices[i].size();
        back_up_file.write(reinterpret_cast<const char*>(&primitive_indices_count), sizeof(uint));
        back_up_file.write(reinterpret_cast<const char*>(m_primitive_indices[i].data()), sizeof(unsigned char) * primitive_indices_count);
    }

    // backup meshlet buffers
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint meshlet_count = (uint)m_meshlets[i].size();
        back_up_file.write(reinterpret_cast<const char*>(&meshlet_count), sizeof(uint));
        back_up_file.write(reinterpret_cast<const char*>(m_meshlets[i].data()), sizeof(S_Meshlet) * meshlet_count);
    }

    back_up_file.close();
    
}

bool Scene::loadSceneFromBackUp(std::string file_path)
{
    std::ifstream back_up_file(file_path, std::ios::binary);
    // early terminate if there is no backup file to load from
    if (!back_up_file.good()) {
        OutputDebugString("FAILED TO LOAD/FIND BACKUP FILE, INITIATE PRE-PROCESSING!\n");
        return false;
    }

    // load scene meta data
    uint scene_object_count;
    back_up_file.read(reinterpret_cast<char*>(&scene_object_count), sizeof(uint));
    back_up_file.read(reinterpret_cast<char*>(&m_mesh_count), sizeof(uint));

    // load scene objects buffer data
    m_scene_objects.resize(scene_object_count);
    back_up_file.read(reinterpret_cast<char*>(m_scene_objects.data()), sizeof(S_SceneObject) * scene_object_count);

    // load vertex buffers
    m_vertices.resize(m_mesh_count);
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint vertex_count; 
        back_up_file.read(reinterpret_cast<char*>(&vertex_count), sizeof(uint));
        m_vertices[i].resize(vertex_count);
        back_up_file.read(reinterpret_cast<char*>(m_vertices[i].data()), sizeof(S_Vertex) * vertex_count);
    }

    // load vertex-indices buffers
    m_vertex_indices.resize(m_mesh_count);
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint vertex_indices_count;
        back_up_file.read(reinterpret_cast<char*>(&vertex_indices_count), sizeof(uint));
        m_vertex_indices[i].resize(vertex_indices_count);
        back_up_file.read(reinterpret_cast<char*>(m_vertex_indices[i].data()), sizeof(uint) * vertex_indices_count);
    }

    // backup morph-indices buffers
    m_morph_indices.resize(m_mesh_count);
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint morph_indices_count;
        back_up_file.read(reinterpret_cast<char*>(&morph_indices_count), sizeof(uint));
        m_morph_indices[i].resize(morph_indices_count);
        back_up_file.read(reinterpret_cast<char*>(m_morph_indices[i].data()), sizeof(uint) * morph_indices_count);
    }

    // backup primitive_indices buffers
    m_primitive_indices.resize(m_mesh_count);
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint primitive_indices_count;
        back_up_file.read(reinterpret_cast<char*>(&primitive_indices_count), sizeof(uint));
        m_primitive_indices[i].resize(primitive_indices_count);
        back_up_file.read(reinterpret_cast<char*>(m_primitive_indices[i].data()), sizeof(unsigned char) * primitive_indices_count);
    }

    // backup meshlet buffers
    m_meshlets.resize(m_mesh_count);
    for (uint i = 0; i < m_mesh_count; i++)
    {
        uint meshlet_count;
        back_up_file.read(reinterpret_cast<char*>(&meshlet_count), sizeof(uint));
        m_meshlets[i].resize(meshlet_count);
        back_up_file.read(reinterpret_cast<char*>(m_meshlets[i].data()), sizeof(S_Meshlet) * meshlet_count);
    }

    back_up_file.close();
    return true;
    
}

void Scene::free()
{
    // reset counter variables
    m_mesh_count = 0;
    m_total_meshlet_count = 0;
    m_totoal_triangle_count = 0;
    m_total_vertex_count = 0;
    m_totalPreProcessingTime = m_totalPreProcessingTime.zero();

    // clear all data vectors
    m_scene_objects.clear();
    for (auto mesh : m_meshes) delete mesh;
    m_meshes.clear();
    m_vertices.clear();
    m_vertex_indices.clear();
    m_primitive_indices.clear();
    m_morph_indices.clear();
    m_meshlets.clear();
}