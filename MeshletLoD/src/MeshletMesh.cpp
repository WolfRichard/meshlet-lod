#include "MeshletMesh.h"
#include "Animator.h"

using namespace DirectX;

MeshletMesh::MeshletMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    parseMesh(assimp_mesh, assimp_scene);

    // LoD
    auto benchmark_time_start = std::chrono::high_resolution_clock::now();
    generateLoD(10);
    auto benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_LoDGenTime = benchmark_time_end - benchmark_time_start;

    // Meshlets
    benchmark_time_start = std::chrono::high_resolution_clock::now();
    generateMeshlets();
    benchmark_time_end = std::chrono::high_resolution_clock::now();
    m_meshletGenTime = benchmark_time_end - benchmark_time_start;
}


void MeshletMesh::parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    m_boundingSphereCentre = float3(0, 0, 0);
    m_boundingSphereRadius = 0;
    // parse raw indices
    std::vector<unsigned int> raw_indices;
    for (uint face = 0; face < assimp_mesh->mNumFaces; face++)
    {
        if (assimp_mesh->mFaces[face].mNumIndices != 3) throw std::exception(); // throw exception if mesh is not triangulated
        for (uint i = 0; i < 3; i++)
        {
            raw_indices.push_back(assimp_mesh->mFaces[face].mIndices[i]);
        }
    }

    // parse raw vertices
    std::vector<CustomVertex> raw_vertices;
    for (uint vert = 0; vert < assimp_mesh->mNumVertices; vert++)
    {
        CustomVertex new_vert;

        new_vert.position = float4(assimp_mesh->mVertices[vert].x,
                                   assimp_mesh->mVertices[vert].y,
                                   assimp_mesh->mVertices[vert].z,
                                   0);

        if (assimp_mesh->HasTextureCoords(0))
        {
            new_vert.uv = float4(assimp_mesh->mTextureCoords[0][vert].x,
                                 assimp_mesh->mTextureCoords[0][vert].y,
                                 0, 0);
        }
        else
        {
            new_vert.uv = float4(0, 0, 0, 0);
        }

        new_vert.normal = float4(assimp_mesh->mNormals[vert].x,
                                 assimp_mesh->mNormals[vert].y,
                                 assimp_mesh->mNormals[vert].z,
                                 0);

        if (assimp_mesh->HasVertexColors(0))
        {
            new_vert.color = float4(assimp_mesh->mColors[0][vert].r,
                                    assimp_mesh->mColors[0][vert].g,
                                    assimp_mesh->mColors[0][vert].b,
                                    assimp_mesh->mColors[0][vert].a);
        }
        else
        {
            new_vert.color = float4(1, 0, 1, 1);
        }

        // set all vertex bone data to default
        for (uint bone = 0; bone < MAX_BONES_PER_VERTEX; bone++) {
            new_vert.bones[bone] = -1;
            new_vert.bone_weights[bone] = 0;
        }
        


        raw_vertices.push_back(new_vert);

        // calc bounding sphere 
        m_boundingSphereCentre.x += new_vert.position.x;
        m_boundingSphereCentre.y += new_vert.position.y;
        m_boundingSphereCentre.z += new_vert.position.z;
    }
    m_boundingSphereCentre.x /= (uint)raw_vertices.size();
    m_boundingSphereCentre.y /= (uint)raw_vertices.size();
    m_boundingSphereCentre.z /= (uint)raw_vertices.size();

    for (uint i = 0; i < (uint)raw_vertices.size(); i++)
    {
        float dx = raw_vertices[i].position.x - m_boundingSphereCentre.x;
        float dy = raw_vertices[i].position.y - m_boundingSphereCentre.y;
        float dz = raw_vertices[i].position.z - m_boundingSphereCentre.z;
        float distSq = dx * dx + dy * dy + dz * dz;
        m_boundingSphereRadius = m_boundingSphereRadius < distSq ? distSq : m_boundingSphereRadius;
    }
    m_boundingSphereRadius = std::sqrt(m_boundingSphereRadius);


    extractAssimpBoneData(raw_vertices, assimp_mesh, assimp_scene);


    // mesh optimize as preperation for lod and meshlet generation
    m_index_count = raw_indices.size();
    std::vector<unsigned int> remap(m_index_count);
    m_vertex_count = meshopt_generateVertexRemap(remap.data(), raw_indices.data(), m_index_count, raw_vertices.data(), raw_vertices.size(), sizeof(CustomVertex));

    m_indices.resize(m_index_count);
    m_vertices.resize(m_vertex_count);

    meshopt_remapIndexBuffer(m_indices.data(), raw_indices.data(), m_index_count, &remap[0]);
    meshopt_remapVertexBuffer(m_vertices.data(), raw_vertices.data(), raw_vertices.size(), sizeof(CustomVertex), &remap[0]);
    meshopt_optimizeVertexCache(m_indices.data(), m_indices.data(), m_index_count, m_vertex_count);

    // check what animations are influencing this mesh
    for (uint a = 0; a < assimp_scene->mNumAnimations; a++)
    {
        bool animationIsCompatibleWithMesh = true;
        if (assimp_mesh->mNumBones >= assimp_scene->mAnimations[a]->mNumChannels)
        {
            for (uint c = 0; c < assimp_scene->mAnimations[a]->mNumChannels; c++)
            {
                std::string boneName = assimp_scene->mAnimations[a]->mChannels[c]->mNodeName.C_Str();
                if (m_BoneInfoMap.find(boneName) == m_BoneInfoMap.end())
                {
                    animationIsCompatibleWithMesh = false;
                }
            }
        }

        if (animationIsCompatibleWithMesh)
        {
            Animator animator;
            Animation animation;
            animation.init(assimp_scene, a, this);
            animator.init(&animation);
            animator.UpdateAnimation(0);

            m_animations.push_back(PreBakedAnimation());
            m_animations.back().boneCount = animator.m_FinalBoneMatrices.size();
            m_animations.back().frameCount = (animation.m_Duration / animation.m_TicksPerSecond) * ANIMATION_FPS;
            m_animations.back().matrices.reserve(m_animations.back().boneCount* m_animations.back().frameCount);
            for (uint f = 0; f < m_animations.back().frameCount; f++)
            {
                m_animations.back().matrices.insert(m_animations.back().matrices.end(), animator.m_FinalBoneMatrices.begin(), animator.m_FinalBoneMatrices.end());
                animator.UpdateAnimation(0.033333333333f);
            }
        }
    }

    
}

void MeshletMesh::generateMeshlets()
{
    const float cone_weight  = 0.0f;
    size_t      max_meshlets = meshopt_buildMeshletsBound(m_indices.size(), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT);
    std::vector<meshopt_Meshlet> m_meshlets(max_meshlets);
    m_meshlet_vertices.resize(max_meshlets * MAX_MESHLET_VERTEX_COUNT);
    m_meshlet_triangles.resize(max_meshlets * MAX_MESHLET_PRIMITIVE_COUNT * 3);

    uint meshlet_count = (uint)meshopt_buildMeshlets(m_meshlets.data(), m_meshlet_vertices.data(), m_meshlet_triangles.data(), m_indices.data(), m_indices.size(), &m_vertices[0].position.x, m_vertices.size(), sizeof(CustomVertex), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT, cone_weight);
    m_meshlets.resize(meshlet_count);

    for (uint i = 0; i < meshlet_count; i++)
    {
        meshopt_Bounds bounds = meshopt_computeMeshletBounds(&m_meshlet_vertices[m_meshlets[i].vertex_offset],
                                                             &m_meshlet_triangles[m_meshlets[i].triangle_offset],
                                                             m_meshlets[i].triangle_count,
                                                             &m_vertices[0].position.x,
                                                             m_vertices.size(),
                                                             sizeof(CustomVertex));


        DrawTask newTask;
        newTask.triangle_count  = m_meshlets[i].triangle_count;
        newTask.triangle_offset = m_meshlets[i].triangle_offset;
        newTask.vertex_count    = m_meshlets[i].vertex_count;
        newTask.vertex_offset   = m_meshlets[i].vertex_offset;


        newTask.culling_info.bounding_sphere_center = float3(bounds.center[0], bounds.center[1], bounds.center[2]);
        newTask.culling_info.bounding_sphere_radius = bounds.radius;
        newTask.culling_info.normal_cone_apex       = float3(bounds.cone_apex[0], bounds.cone_apex[1], bounds.cone_apex[2]);
        newTask.culling_info.normal_cone_axis       = float3(bounds.cone_axis[0], bounds.cone_axis[1], bounds.cone_axis[2]);
        newTask.culling_info.normal_cone_cutoff     = bounds.cone_cutoff;
        newTask.culling_info.byte_alignement        = 0.f;


        m_draw_tasks.push_back(newTask);
    }
}

void MeshletMesh::generateLoD(uint number_of_levels) 
{
    for (uint resolution_level = 1; resolution_level < number_of_levels; resolution_level++)
    {
        size_t target_index_count = static_cast<size_t>((uint)std::max(m_indices.size() / (uint)std::pow(2, resolution_level), 3ull));
        std::vector<uint> simplified_indices(m_indices.size());

        // Simplify the mesh
        float lod_error = 0.f;
        size_t simplified_count = meshopt_simplifySloppy(
            simplified_indices.data(), m_indices.data(), m_indices.size(),
            reinterpret_cast<const float*>(m_vertices.data()), m_vertices.size(), sizeof(CustomVertex),
            target_index_count, FLT_MAX, &lod_error
        );

        // Resize the index vector to the actual simplified size
        simplified_indices.resize(simplified_count);
        //m_indices = simplified_indices;
    }
}

void MeshletMesh::addVertexBoneData(CustomVertex& vertex, int boneID, float weight)
{
    for (uint bone = 0; bone < MAX_BONES_PER_VERTEX; bone++) 
    {
        if (vertex.bones[bone] < 0) // check if bone slot is free
        {
            vertex.bones[bone] = boneID;
            vertex.bone_weights[bone] = weight;
            return;
        }
    }
}

// TODO: CHech if there is a solution without a map is possible by just using a vector to store the bones and use their index position as boneID (just a vector of their matrices)

void MeshletMesh::extractAssimpBoneData(std::vector<CustomVertex>& vertices, aiMesh* mesh, const aiScene* scene)
{
    for (uint boneIndex = 0; boneIndex < mesh->mNumBones; ++boneIndex)
    {
        int boneID = -1;
        std::string boneName = mesh->mBones[boneIndex]->mName.C_Str();
        if (m_BoneInfoMap.find(boneName) == m_BoneInfoMap.end())
        {
            BoneInfo newBoneInfo;
            newBoneInfo.id = m_BoneCounter;
            newBoneInfo.transform_matrix = assimpToFloat4x4Matrix(mesh->mBones[boneIndex]->mOffsetMatrix);
            m_BoneInfoMap[boneName] = newBoneInfo;
            boneID = m_BoneCounter;
            m_BoneCounter++;
        }
        else
        {
            boneID = m_BoneInfoMap[boneName].id;
        }
        assert(boneID != -1);
        auto weights = mesh->mBones[boneIndex]->mWeights;
        int numWeights = mesh->mBones[boneIndex]->mNumWeights;

        for (int weightIndex = 0; weightIndex < numWeights; ++weightIndex)
        {
            int vertexId = weights[weightIndex].mVertexId;
            float weight = weights[weightIndex].mWeight;
            assert(vertexId <= vertices.size());
            addVertexBoneData(vertices[vertexId], boneID, weight);
        }
    }
}
