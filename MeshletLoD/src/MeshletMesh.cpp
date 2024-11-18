#include "MeshletMesh.h"

MeshletMesh::MeshletMesh(aiMesh* assimp_mesh)
{
    parseMesh(assimp_mesh);
    generateMeshlets();
}


void MeshletMesh::parseMesh(aiMesh* assimp_mesh)
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


    // mesh optimize as preperation for lod and meshlet generation
    m_index_count = raw_indices.size();
    std::vector<unsigned int> remap(m_index_count);
    m_vertex_count = meshopt_generateVertexRemap(remap.data(), raw_indices.data(), m_index_count, raw_vertices.data(), raw_vertices.size(), sizeof(CustomVertex));

    m_indices.resize(m_index_count);
    m_vertices.resize(m_vertex_count);

    meshopt_remapIndexBuffer(m_indices.data(), raw_indices.data(), m_index_count, &remap[0]);
    meshopt_remapVertexBuffer(m_vertices.data(), raw_vertices.data(), raw_vertices.size(), sizeof(CustomVertex), &remap[0]);
    meshopt_optimizeVertexCache(m_indices.data(), m_indices.data(), m_index_count, m_vertex_count);
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

