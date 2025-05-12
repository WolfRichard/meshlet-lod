#include "Mesh.h"

#include <deque>
#include <metis.h>
#include <unordered_set>

#include <random>

#include "Windows.h"


using namespace DirectX;

Mesh::Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    OutputDebugString("Started Mesh Initialisation\n");
    parseMesh(assimp_mesh, assimp_scene);
    OutputDebugString("Done parsing\n");
    generateLeafMeshlets();
    OutputDebugString("Generated Leaf Meshlets\n");
    buildMeshletHierachy();
    OutputDebugString("Finished Hierarchy Generation\n");
}


void Mesh::parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    // initialize bounding sphere approximation with zero
    // --> set bounding sphere radius ontop of average vertex position
    // --> set radius = furthest vertex distance from center
    m_bounding_sphere.center = float3(0, 0, 0);
    m_bounding_sphere.radius = 0;

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
    std::vector<S_Vertex> raw_vertices;
    for (uint vert = 0; vert < assimp_mesh->mNumVertices; vert++)
    {
        S_Vertex new_vert;

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

        // calculate average vertex position 
        m_bounding_sphere.center.x += new_vert.position.x;
        m_bounding_sphere.center.y += new_vert.position.y;
        m_bounding_sphere.center.z += new_vert.position.z;
    }
    m_bounding_sphere.center.x /= (uint)raw_vertices.size();
    m_bounding_sphere.center.y /= (uint)raw_vertices.size();
    m_bounding_sphere.center.z /= (uint)raw_vertices.size();

    // calculate conservative bounding sphere radius
    for (uint i = 0; i < (uint)raw_vertices.size(); i++)
    {
        float dx = raw_vertices[i].position.x - m_bounding_sphere.center.x;
        float dy = raw_vertices[i].position.y - m_bounding_sphere.center.y;
        float dz = raw_vertices[i].position.z - m_bounding_sphere.center.z;
        float distSq = dx * dx + dy * dy + dz * dz;
        m_bounding_sphere.radius = m_bounding_sphere.radius < distSq ? distSq : m_bounding_sphere.radius;
    }
    m_bounding_sphere.radius = std::sqrt(m_bounding_sphere.radius);

    // mesh optimize as preperation for lod and meshlet generation
    std::vector<unsigned int> remap(raw_indices.size());
    size_t vertex_count = meshopt_generateVertexRemap(remap.data(), raw_indices.data(), raw_indices.size(), raw_vertices.data(), raw_vertices.size(), sizeof(S_Vertex));

    m_original_indices.resize(raw_indices.size());
    m_vertices.resize(vertex_count);

    meshopt_remapIndexBuffer(m_original_indices.data(), raw_indices.data(), raw_indices.size(), &remap[0]);
    meshopt_remapVertexBuffer(m_vertices.data(), raw_vertices.data(), raw_vertices.size(), sizeof(S_Vertex), &remap[0]);
    meshopt_optimizeVertexCache(m_original_indices.data(), m_original_indices.data(), raw_indices.size(), vertex_count);

    // simplify mesh if triangle count is larger than maximum to prevent excessive pre-compute times
    if (m_maximum_mesh_triangle_count)
    {
        size_t target_index_count = min(m_maximum_mesh_triangle_count * 3, m_original_indices.size());
        std::vector<uint> simplified_indices(m_original_indices);
        float lod_error = 0.0f;
        size_t simplifiedIndexCount = meshopt_simplify(
            simplified_indices.data(), m_original_indices.data(), m_original_indices.size(),
            reinterpret_cast<const float*>(&m_vertices.data()[0].position.x), m_vertices.size(), sizeof(S_Vertex),
            target_index_count, FLT_MAX, meshopt_SimplifyLockBorder, &lod_error); // & meshopt_SimplifyErrorAbsolute bitmask option leads to cracks??? why??
        simplified_indices.resize(simplifiedIndexCount);
        m_original_indices = simplified_indices;
    }
}

void Mesh::generateLeafMeshlets()
{
    // generate meshlets using the meshoptimizer library
    const float cone_weight = 0.0f;
    size_t      max_meshlets = meshopt_buildMeshletsBound(m_original_indices.size(), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT);
    std::vector<meshopt_Meshlet> meshlets(max_meshlets);
    m_vertex_indices.resize(max_meshlets * MAX_MESHLET_VERTEX_COUNT);
    m_primitive_indices.resize(max_meshlets * MAX_MESHLET_PRIMITIVE_COUNT * 3);
    uint meshlet_count = (uint)meshopt_buildMeshlets(meshlets.data(), m_vertex_indices.data(), m_primitive_indices.data(), m_original_indices.data(), m_original_indices.size(), &m_vertices[0].position.x, m_vertices.size(), sizeof(S_Vertex), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT, cone_weight);
    meshlets.resize(meshlet_count);

    // set the first meshlet to be an empty meshlet
    // this will be used to balance holes in groups that cant be filled up to 4 meshlets, always at id 0
    S_Meshlet dummyMeshlet; 
    dummyMeshlet.triangle_count = 0;
    dummyMeshlet.triangle_offset = 0;
    dummyMeshlet.vertex_count = 0;
    dummyMeshlet.vertex_offset = 0;
    dummyMeshlet.base_error = 0;
    dummyMeshlet.discrete_level_of_detail = 0;
    dummyMeshlet.simplification_error = FLT_MAX;
    dummyMeshlet.bounding_sphere = { float3(0, 0, 0), 0 };
    for (uint i = 0; i < GROUP_MERGE_COUNT; i++)
        dummyMeshlet.child_meshlets[i] = 0;
    dummyMeshlet.child_count = 0;

    m_meshlets.push_back(dummyMeshlet);

    // parse meshoptimizer meshlets into own data structures
    for (uint i = 0; i < meshlet_count; i++)
    {
        S_Meshlet newMeshlet;
        newMeshlet.triangle_count = meshlets[i].triangle_count;
        newMeshlet.triangle_offset = meshlets[i].triangle_offset;
        newMeshlet.vertex_count = meshlets[i].vertex_count;
        newMeshlet.vertex_offset = meshlets[i].vertex_offset;
        newMeshlet.base_error = 0;
        newMeshlet.discrete_level_of_detail = 0;
        newMeshlet.simplification_error = FLT_MAX;
        newMeshlet.simplified_group_bounds = { float3(0, 0, 0), 0 };
        for (uint i = 0; i < GROUP_MERGE_COUNT; i++)
            newMeshlet.child_meshlets[i] = 0;
        newMeshlet.child_count = 0;

        // generate and parse bounding data for the meshlets
        meshopt_Bounds bounds = meshopt_computeMeshletBounds(&m_vertex_indices[newMeshlet.vertex_offset],
            &m_primitive_indices[newMeshlet.triangle_offset],
            newMeshlet.triangle_count,
            &m_vertices[0].position.x,
            m_vertices.size(),
            sizeof(S_Vertex)
        );
        newMeshlet.bounding_sphere.center = float3(bounds.center[0], bounds.center[1], bounds.center[2]);
        newMeshlet.bounding_sphere.radius = bounds.radius;

        // set newly generated meshlets as the currently processed level of the hierarchy, so that the other functions know what to work on
        m_current_hierarchy_top_level_meshlets.push_back((int)m_meshlets.size());
        m_meshlets.push_back(newMeshlet);
    }

    // initialize morph indices to point at themselfes (no morphing) when first creating base meshlets
    m_morph_indices.resize(m_vertex_indices.size(), 0);
    for (uint m = 0; m < m_current_hierarchy_top_level_meshlets.size(); m++)
    {
        S_Meshlet& current_simplified_meshlet = m_meshlets[m_current_hierarchy_top_level_meshlets[m]];
        for (unsigned char vi = 0; vi < current_simplified_meshlet.vertex_count; vi++)
        {
            m_morph_indices[current_simplified_meshlet.vertex_offset + vi] = current_simplified_meshlet.vertex_offset + vi;
        }
    }
}

void Mesh::buildMeshletHierachy() 
{
    // group and simplify as long as there are enough meshlets to do so
    while ((uint)m_current_hierarchy_top_level_meshlets.size() >= GROUP_MERGE_COUNT * 2)
    {
        groupMeshlets();
        simplifiyTopLevelGroups();
    }
   
    assert(m_current_hierarchy_top_level_groups.size() == 2);
    OutputDebugString(("\nNumber of Meshlets: " + std::to_string(m_current_hierarchy_top_level_meshlets.size()) + "\n").c_str());
    OutputDebugString(("Number of Groups: " + std::to_string(1) + "\n").c_str());
    
    // edge case handeling for root nodes
    finalTopLevelMeshletGrouping();
    simplifiyTopLevelGroups();

    // store root meshlet id's inside the dummy meshlet (maybe not the most elegant solution but it saves the need for an extra buffer)
    m_meshlets[0].child_meshlets[0] = (uint)m_meshlets.size() - 2;
    m_meshlets[0].child_meshlets[1] = (uint)m_meshlets.size() - 1;

    //clear uneeded data
    m_current_hierarchy_top_level_groups.clear();
    m_current_hierarchy_top_level_meshlets.clear();

    // print debug info
    OutputDebugString(("\nTotal number of Meshlets: " + std::to_string(m_meshlets.size() - 1) + "\n").c_str());
    OutputDebugString(("Total number of Groups: " + std::to_string(m_meshlet_groups.size()) + "\n").c_str());
    OutputDebugString(("Hierarchy Tree Depth: " + std::to_string(m_hierarchy_per_level_group_count.size()) + "\n").c_str());
}



void Mesh::groupMeshlets()
{
    OutputDebugString(("\nNumber of Meshlets: " + std::to_string(m_current_hierarchy_top_level_meshlets.size()) + "\n").c_str());

    // build METIS graph representation of meshlet connectivity to each other
    idx_t MATIS_numNodes = (int)m_current_hierarchy_top_level_meshlets.size();
    std::vector<idx_t> MATIS_vertexWeights(MATIS_numNodes, 100);
    idx_t MATIS_numEdges = 0;
    std::vector<idx_t> MATIS_nodeAdjacencyOffsets;
    std::vector<idx_t> MATIS_adjacencyList;
    std::vector<idx_t> MATIS_edgeWeights; // represents the number of shared edges between each meshlet

    // extract all boundary edges from the current top level meshlets
    std::vector<std::unordered_set<std::pair<uint, uint>, PairHash>> meshlet_edges;
    for (auto meshlet_index : m_current_hierarchy_top_level_meshlets)
    {
        meshlet_edges.push_back(extractBoundaryEdgeSet(m_meshlets[meshlet_index]));
    }
    std::vector<std::vector<uint>> connectivity_matrix;

    // count shared edges between every meshlet
    for (uint a = 0; a < m_current_hierarchy_top_level_meshlets.size(); a++)
    {
        uint meshlet_A_index = m_current_hierarchy_top_level_meshlets[a];
        connectivity_matrix.push_back(std::vector<uint>());
        MATIS_nodeAdjacencyOffsets.push_back((int)MATIS_adjacencyList.size());
        for (uint b = 0; b < m_current_hierarchy_top_level_meshlets.size(); b++)
        {
            uint meshlet_B_index = m_current_hierarchy_top_level_meshlets[b];
            if (meshlet_A_index == meshlet_B_index)
            {
                connectivity_matrix.back().push_back(0);
            }
            else 
            { 
                uint shared_edge_count = count_shared_edges(meshlet_edges[a], meshlet_edges[b]);
                connectivity_matrix.back().push_back(shared_edge_count);
                if (shared_edge_count)
                {
                    MATIS_adjacencyList.push_back(b);
                    MATIS_edgeWeights.push_back(shared_edge_count);
                    if (a > b) // prevents edges from beeing counted double (only counts one way)
                        MATIS_numEdges++;
                }
                    
            }
        }
    }
    MATIS_nodeAdjacencyOffsets.push_back((int)MATIS_adjacencyList.size());
    idx_t MATIS_numPartitions = (uint)((m_current_hierarchy_top_level_meshlets.size() / GROUP_MERGE_COUNT));// + GROUP_MERGE_COUNT - 1) / GROUP_MERGE_COUNT);
    std::vector<idx_t> MATIS_outputPartitionsArray(MATIS_numNodes);
    idx_t METIS_options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(METIS_options); // Use default options
    METIS_options[METIS_OPTION_UFACTOR] = 100;  // Strictest balance (but may increase edge cut)
    idx_t MATIS_ncon = 1;
    idx_t MATIS_edgeCut;
    int result = METIS_PartGraphKway(&MATIS_numNodes,                           // number of nodes in the graph (= numbers of currently processed meshlets)
                                          &MATIS_ncon,                               // number of balancing constraints (not needed thus default = 1)
                                          MATIS_nodeAdjacencyOffsets.data(),         // offsets for each nodes neighbor list
                                          MATIS_adjacencyList.data(),                // list of edges of the graph represented as array of neighbors of each node (so basicly each edge is stored twice)
                                          MATIS_vertexWeights.data(),                                   // vertex weights (not needed thus disabled = nullptr)
                                          nullptr,                                   // vertex sizes (also not necessary = nullptr)
                                          MATIS_edgeWeights.data(),                  // edge weights list (weight for each neighbor connection of the adjacencyList, so gain every weight s stored twice)
                                          &MATIS_numPartitions,                      // number of partitions that the graph should be devided in (how many groups of 4 meshlets we want)
                                          nullptr,                                   // desired weight distribution between partitions (we want equal distribution so default = nullptr)
                                          nullptr,                                   // balance tolerance for each constraint (not needed thus = nullptr)
                                          METIS_options,                             // additional options for graph partitioning (we stick to default settings)
                                          &MATIS_edgeCut,                            // returns the number of edges that the resulting partitoning has cut to divide the graph
                                          MATIS_outputPartitionsArray.data()         // returns the group id for each node that the meshlet is part of
                                          );

    std::vector<uint> groupMeshletCounts(MATIS_numPartitions); // holds the number of meshlets that are part of each group
    for (int k = 0; k < MATIS_numPartitions; k++)
    {
        uint meshletCounter = 0;
        for (auto assignedGroupIndex : MATIS_outputPartitionsArray)
            if (k == assignedGroupIndex)
                meshletCounter++;
        groupMeshletCounts[k] = meshletCounter;
    }

    OutputDebugString(("Number of Groups: " + std::to_string(MATIS_numPartitions) + "\n").c_str());

    // mark newly gnerated meshlet grous as current top level 
    m_current_hierarchy_top_level_groups.clear();
    for (int k = 0; k < MATIS_numPartitions; k++)
    {
        m_current_hierarchy_top_level_groups.push_back((uint)m_meshlet_groups.size());
        m_meshlet_groups.push_back(getDefaultMeshletGroup());
        S_MeshletGroup& current_group = m_meshlet_groups.back();
        current_group.hierarchy_tree_depth = (uint)m_hierarchy_per_level_group_count.size();

        // set group child meshlets
        current_group.meshlet_count = 0;
        for (int j = 0; j < MATIS_numNodes; j++)
        {
            if (k == MATIS_outputPartitionsArray[j])
            {
                assert(current_group.meshlet_count < (sizeof(current_group.meshlets) / 4));
                current_group.meshlets[current_group.meshlet_count++] = m_current_hierarchy_top_level_meshlets[j];
            }
        }
    }
    m_hierarchy_per_level_group_count.push_back((uint)m_current_hierarchy_top_level_groups.size());
}


std::pair<uint, uint> Mesh::sortEdgeIndices(uint v0, uint v1) 
{
    // always push lesser index to be first
    return (v0 < v1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
}


std::unordered_set<std::pair<uint, uint>, PairHash> Mesh::extractEdges(S_Meshlet meshlet) 
{
    std::unordered_set<std::pair<uint, uint>, PairHash> edge_set;

    // extract all 3 edges for every triangle and return the reult as an unordered_set for faster lookup
    for (uint i = 0; i < meshlet.triangle_count; i++) 
    {
        uint vertex_index_0 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 0]];
        uint vertex_index_1 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 1]];
        uint vertex_index_2 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 2]];

        edge_set.insert(sortEdgeIndices(vertex_index_0, vertex_index_1));
        edge_set.insert(sortEdgeIndices(vertex_index_1, vertex_index_2));
        edge_set.insert(sortEdgeIndices(vertex_index_2, vertex_index_0));
    }

    return edge_set;
}


std::unordered_set<std::pair<uint, uint>, PairHash> Mesh::extractBoundaryEdgeSet(S_Meshlet meshlet)
{
    std::unordered_set<std::pair<uint, uint>, PairHash> edge_set;
    std::unordered_map<std::pair<uint, uint>, int, PairHash> edgeCount;

    // count the occurance of every edge
    for (uint i = 0; i < meshlet.triangle_count; i++)
    {
        uint vertex_index_0 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 0]];
        uint vertex_index_1 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 1]];
        uint vertex_index_2 = m_vertex_indices[meshlet.vertex_offset + m_primitive_indices[meshlet.triangle_offset + i * 3 + 2]];

        edgeCount[sortEdgeIndices(vertex_index_0, vertex_index_1)]++;
        edgeCount[sortEdgeIndices(vertex_index_1, vertex_index_2)]++;
        edgeCount[sortEdgeIndices(vertex_index_2, vertex_index_0)]++;
    }

    // boundary edges = edges that only occur once
    for (const auto& [edge, count] : edgeCount) {
        if (count == 1) {
            edge_set.insert(edge);
        }
    }

    return edge_set;
}


uint Mesh::count_shared_edges(const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesA, const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesB) 
{
    // go over every edge of set A and try to find it in set B
    int shared_count = 0;
    for (const auto& edge : edgesA) 
    {
        if (edgesB.find(edge) != edgesB.end()) 
        {
            shared_count++;
        }
    }
    return shared_count;
}


std::vector<std::pair<uint, uint>> Mesh::extractBoundaryEdges(std::vector<uint>& indices)
{
    std::unordered_map<std::pair<uint, uint>, int, PairHash> edgeCount;

    // count occurrences of each edge
    for (uint tri = 0; tri < (uint)indices.size() / 3; tri++)
    {
        edgeCount[sortEdgeIndices(indices[tri * 3 + 0], indices[tri * 3 + 1])]++;
        edgeCount[sortEdgeIndices(indices[tri * 3 + 1], indices[tri * 3 + 2])]++;
        edgeCount[sortEdgeIndices(indices[tri * 3 + 2], indices[tri * 3 + 0])]++;
    }

    // boundary edges = edges that only occur once
    std::vector<std::pair<uint, uint>> boundaryEdges;
    for (const auto& [edge, count] : edgeCount) {
        if (count == 1) { 
            boundaryEdges.push_back(edge);
        }
    }

    return boundaryEdges;
}


void Mesh::simplifiyTopLevelGroups()
{
    // clear the top level of meshlets as it will be later replaced by their simplified version
    m_current_hierarchy_top_level_meshlets.clear(); 

    // process each group individually
    for (uint g = 0; g < m_current_hierarchy_top_level_groups.size(); g++) 
    {
        OutputDebugString(("\t -> simplifiying group[" + std::to_string(g) + "]\n").c_str());
        
        S_MeshletGroup& current_group = m_meshlet_groups[m_current_hierarchy_top_level_groups[g]];

        // combined triangles of all meshlets (in reference to the original vertex buffer)
        std::vector<uint> merged_deduplicated_indices; 
        std::vector<S_BoundingSphere> meshlet_bounding_spheres;
        for (uint m = 0; m < current_group.meshlet_count; m++)
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
            current_meshlet.group_id = m_current_hierarchy_top_level_groups[g];
            
            current_group.bounding_sphere.center.x += current_meshlet.bounding_sphere.center.x;
            current_group.bounding_sphere.center.y += current_meshlet.bounding_sphere.center.y;
            current_group.bounding_sphere.center.z += current_meshlet.bounding_sphere.center.z;
            
            meshlet_bounding_spheres.push_back(current_meshlet.bounding_sphere);

            for (uint i = 0; i < current_meshlet.triangle_count * 3; i++) 
            {
                merged_deduplicated_indices.push_back(m_vertex_indices[current_meshlet.vertex_offset + m_primitive_indices[current_meshlet.triangle_offset + i]]);
            }
        }

        // average out the bounding sphere centers of chid meshlets to get approximateion for group center
        current_group.bounding_sphere.center.x /= current_group.meshlet_count;
        current_group.bounding_sphere.center.y /= current_group.meshlet_count;
        current_group.bounding_sphere.center.z /= current_group.meshlet_count;

        // compute conservative bounding sphere radius by looking for furthest outgoing inner meshlet bounding sphere
        current_group.bounding_sphere.radius = 0;
        for (uint m = 0; m < current_group.meshlet_count; m++)
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
            float dx = current_group.bounding_sphere.center.x - current_meshlet.bounding_sphere.center.x;
            float dy = current_group.bounding_sphere.center.y - current_meshlet.bounding_sphere.center.y;
            float dz = current_group.bounding_sphere.center.z - current_meshlet.bounding_sphere.center.z;

            float dist = sqrtf(dx * dx + dy * dy + dz * dz) + current_meshlet.bounding_sphere.radius;
            current_group.bounding_sphere.radius = current_group.bounding_sphere.radius < dist ? dist : current_group.bounding_sphere.radius;
        }
         
        // optional alternative more accurate boundingsphere calculation (vertex based instead of based on child bounds)
        //computeGroupBoundingSphere(current_group)


        // simplify the meshlet group to have enough space to hold 2 meshlets
        std::vector<uint> simplified_indices(merged_deduplicated_indices.size());
        float lod_error = 0.0f;       

        // stores what vertex was simplified into what other vertex
        std::vector<uint> vertex_remap(m_vertices.size(), -1);
        for (int i = 0; i < m_vertices.size(); i++)
        {
            vertex_remap[i] = i;
        }
        std::vector<int> per_vertex_morph_index_remap(m_vertices.size(), -1);

        // simplify group geometry but keep edge vertices locked
        size_t target_index_count = GROUP_SPLIT_COUNT * MAX_MESHLET_VERTEX_COUNT * 3; // TODO: should be "GROUP_SPLIT_COUNT * MAX_MESHLET_PRIMITIVE_COUNT * 3", but that would blow past the vertex limit and generate to many meshlets
        size_t simplifiedIndexCount = meshopt_simplify_tracking(
            simplified_indices.data(), merged_deduplicated_indices.data(), merged_deduplicated_indices.size(),
            reinterpret_cast<const float*>(&m_vertices.data()[0].position.x), m_vertices.size(), sizeof(S_Vertex),
            target_index_count, FLT_MAX, meshopt_SimplifyLockBorder /* | meshopt_SimplifySparse */ | meshopt_SimplifyErrorAbsolute, &lod_error, &vertex_remap);
        simplified_indices.resize(simplifiedIndexCount);
        
        // set base error line for simplified meshlets aas highest error from base meshlets
        float collective_base_meshlet_error = 0;
        for (uint m = 0; m < current_group.meshlet_count; m++)
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
            collective_base_meshlet_error = max(current_meshlet.base_error, collective_base_meshlet_error);
        }
        collective_base_meshlet_error += lod_error;

        // set simplified error for base meshlets
        for (uint m = 0; m < current_group.meshlet_count; m++)
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
            current_meshlet.simplification_error = collective_base_meshlet_error;
            current_meshlet.simplified_group_bounds = current_group.bounding_sphere;
        }

        // generate new meshlets from simplified geometry
        const float cone_weight = 0.0f;
        size_t      max_meshlets = meshopt_buildMeshletsBound(simplified_indices.size(), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT);
        std::vector<meshopt_Meshlet> meshlets(max_meshlets);
        std::vector<uint> vertex_indices(max_meshlets * MAX_MESHLET_VERTEX_COUNT);
        std::vector<unsigned char> primitive_indices(max_meshlets * MAX_MESHLET_PRIMITIVE_COUNT * 3);
        uint meshlet_count = (uint)meshopt_buildMeshlets(meshlets.data(), vertex_indices.data(), primitive_indices.data(), simplified_indices.data(), simplified_indices.size(), &m_vertices[0].position.x, m_vertices.size(), sizeof(S_Vertex), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT, cone_weight);
        meshlets.resize(meshlet_count);
        assert(meshlet_count == 2);

        // add new vertex/primitive-indicesbuffers to the main buffer and offset meshlet meta data accordingly
        uint previous_primitive_indices_size = (uint)m_primitive_indices.size();
        uint previous_vertex_indices_size    = (uint)m_vertex_indices.size();

        m_vertex_indices.insert(m_vertex_indices.end(), vertex_indices.begin(), vertex_indices.end());
        m_primitive_indices.insert(m_primitive_indices.end(), primitive_indices.begin(), primitive_indices.end());

        // parse each newly generated meshlet into own data format
        for (uint i = 0; i < meshlet_count; i++)
        {
            S_Meshlet newMeshlet;
            newMeshlet.triangle_count = meshlets[i].triangle_count;
            newMeshlet.triangle_offset = meshlets[i].triangle_offset + previous_primitive_indices_size;
            newMeshlet.vertex_count = meshlets[i].vertex_count;
            newMeshlet.vertex_offset = meshlets[i].vertex_offset + previous_vertex_indices_size;
            newMeshlet.base_error = collective_base_meshlet_error;
            newMeshlet.discrete_level_of_detail = (uint)m_hierarchy_per_level_group_count.size();
            newMeshlet.simplification_error = FLT_MAX;
            newMeshlet.simplified_group_bounds = { float3(0, 0, 0), 0 };
            for (uint i = 0; i < GROUP_MERGE_COUNT; i++)
                newMeshlet.child_meshlets[i] = 0;
            newMeshlet.child_count = 0;

            meshopt_Bounds bounds = meshopt_computeMeshletBounds(&m_vertex_indices[newMeshlet.vertex_offset],
                &m_primitive_indices[newMeshlet.triangle_offset],
                newMeshlet.triangle_count,
                &m_vertices[0].position.x,
                m_vertices.size(),
                sizeof(S_Vertex)
            );
            newMeshlet.bounding_sphere.center = float3(bounds.center[0], bounds.center[1], bounds.center[2]);
            newMeshlet.bounding_sphere.radius = bounds.radius;
            newMeshlet.group_id = 0;

            // bounding sphere has to be unique for all meshlets of a group, so that all meshlets perform the same tree cut decision
            newMeshlet.bounding_sphere = current_group.bounding_sphere; 

            // store parsed meshlet in the buffer and push its id into the top level vector
            m_current_hierarchy_top_level_meshlets.push_back((uint)m_meshlets.size());
            current_group.simplified_meshlets[i] = (uint)m_meshlets.size();
            m_meshlets.push_back(newMeshlet);
        }

        // distribute child connections evenly between simplyfied meshlet parrents to balance tree out
        for (uint bm = 0; bm < current_group.meshlet_count; bm++)
        {
            S_Meshlet& simpl_meshlet = m_meshlets[current_group.simplified_meshlets[bm % GROUP_SPLIT_COUNT]];
            simpl_meshlet.child_meshlets[simpl_meshlet.child_count] = current_group.meshlets[bm];
            simpl_meshlet.child_count++;
        }

        // transform morph_target_index_indices  destination (value) entries to reference into the m_vertex_indice buffer instead
        for (uint m = 0; m < GROUP_SPLIT_COUNT; m++)
        {
            S_Meshlet& current_simplified_meshlet = m_meshlets[current_group.simplified_meshlets[m]];
            for (unsigned char svi = 0; svi < current_simplified_meshlet.vertex_count; svi++)
            {
                uint simplified_vertex_index = m_vertex_indices[current_simplified_meshlet.vertex_offset + svi];

                for (uint i = 0; i < m_vertices.size(); i++)
                {
                    if (vertex_remap[i] == simplified_vertex_index)
                        per_vertex_morph_index_remap[i] = current_simplified_meshlet.vertex_offset + svi;
                }
            }
        }
        m_morph_indices.resize(m_vertex_indices.size(), 0); // resize in case new meshlets got added in previous hierarchy expansion
        for (uint m = 0; m < current_group.meshlet_count; m++) // go over every meshlet of the current group
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];

            for (uint vi = 0; vi < current_meshlet.vertex_count; vi++)
            {
                uint vertex_index = m_vertex_indices[current_meshlet.vertex_offset + vi];
                if (per_vertex_morph_index_remap[vertex_index] >= 0)
                {
                    m_morph_indices[current_meshlet.vertex_offset + vi] = per_vertex_morph_index_remap[vertex_index];
                }
            }
        }
        
        for (uint m = 0; m < GROUP_SPLIT_COUNT; m++) // go over every simplified meshlet of the current group
        {
            S_Meshlet& current_simplified_meshlet = m_meshlets[current_group.simplified_meshlets[m]];
            for (unsigned char vi = 0; vi < current_simplified_meshlet.vertex_count; vi++) 
            {
                m_morph_indices[current_simplified_meshlet.vertex_offset + vi] = current_simplified_meshlet.vertex_offset + vi;
            }
        }
    }
}

void Mesh::finalTopLevelMeshletGrouping()
{
    // handle edge case of hierarhy root note group
    m_current_hierarchy_top_level_groups.clear();
    m_hierarchy_root_group = (uint)m_meshlet_groups.size();
    m_current_hierarchy_top_level_groups.push_back((uint)m_meshlet_groups.size());
    m_meshlet_groups.push_back(getDefaultMeshletGroup());
    S_MeshletGroup& current_group = m_meshlet_groups.back();
    current_group.hierarchy_tree_depth = (uint)m_hierarchy_per_level_group_count.size();

    // set group child meshlets
    current_group.meshlet_count = 0;
    for (uint m : m_current_hierarchy_top_level_meshlets)
    {
        assert(current_group.meshlet_count < (sizeof(current_group.meshlets) / 4));
        current_group.meshlets[current_group.meshlet_count++] = m;
    }
    m_hierarchy_per_level_group_count.push_back((uint)m_current_hierarchy_top_level_groups.size());
}

S_MeshletGroup Mesh::getDefaultMeshletGroup()
{
    // initialize every struct variable with either 0 or a conservative value that would result in an empty group
    S_MeshletGroup ret;
    ret.bounding_sphere.center = float3(0, 0, 0);
    ret.bounding_sphere.radius = 0;
    ret.childCount = 0;
    ret.meshlet_count = 0;
    for (uint i = 0; i < sizeof(ret.children) / 4; i++)
    {
        ret.children[i] = -1;
    }
    for (uint j = 0; j < sizeof(ret.meshlets) / 4; j++)
    {
        ret.meshlets[j] = 0;
    }
    for (uint k = 0; k < sizeof(ret.simplified_meshlets) / 4; k++)
    {
        ret.simplified_meshlets[k] = 0;
    }
    ret.parent = -1;
    ret.hierarchy_tree_depth = 0;
    return ret;
}

void Mesh::printTreeFromTop(uint current_group_index) 
{
    // recursive tree traversal for debug printing group hierarchy 
    S_MeshletGroup& current_group = m_meshlet_groups[current_group_index];
    OutputDebugString(("\nGroup[" + std::to_string(current_group_index) + "] --> Children: ").c_str());
    for (uint i = 0; i < current_group.childCount; i++)
        OutputDebugString(("[" + std::to_string(current_group.children[i]) + "], ").c_str());
    for (uint i = 0; i < current_group.childCount; i++)
        printTreeFromTop(current_group.children[i]);
}

void Mesh::printAllMeshlets()
{
    // itterate over every meshlet and debug print its stats
    for (uint m = 0; m < m_meshlets.size(); m++)
    {
        S_Meshlet& current_meshlet = m_meshlets[m];
        OutputDebugString(("\nMeshlet[" + std::to_string(m) + "] - vertex_count: " + std::to_string(current_meshlet.vertex_count) + ", vertex_offset: "
            + std::to_string(current_meshlet.vertex_offset) + ", triangle_count: " + std::to_string(current_meshlet.triangle_count) + " triangle_offset: " 
            + std::to_string(current_meshlet.triangle_offset)).c_str());
    }
}

S_BoundingSphere Mesh::computeGroupBoundingSphere(S_MeshletGroup& meshlet_group)
{
    S_BoundingSphere result;
    result.center = float3(0, 0, 0);
    result.radius = 0;
    if (!(meshlet_group.childCount)) return result;

    // estimate center by averaging out vertex positions
    uint counter = 0;
    for (uint m = 0; m < meshlet_group.meshlet_count; m++)
    {
        S_Meshlet& current_meshlet = m_meshlets[meshlet_group.meshlets[m]];
        for (uint i = 0; i < current_meshlet.vertex_count; i++)
        {
            result.center.x += m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.x;
            result.center.y += m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.y;
            result.center.z += m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.z;
            counter++;
        }
    }
    result.center.x /= counter;
    result.center.y /= counter;
    result.center.z /= counter;


    // calculate conservative radius from estimated center
    float max_dist2 = 0;
    for (uint m = 0; m < meshlet_group.meshlet_count; m++)
    {
        S_Meshlet& current_meshlet = m_meshlets[meshlet_group.meshlets[m]];
        for (uint i = 0; i < current_meshlet.vertex_count; i++)
        {
            float dx = result.center.x - m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.x;
            float dy = result.center.x - m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.y;
            float dz = result.center.x - m_vertices[m_vertex_indices[current_meshlet.vertex_offset + i]].position.z;
            float dist2 = dx * dx + dy * dy + dz * dz;

            max_dist2 = max(max_dist2, dist2);
        }
    }
    result.radius = sqrtf(max_dist2);

    return result;
}


int Mesh::randomInt(int min, int max) {
    static std::random_device rd;  // Seed
    static std::mt19937 gen(rd()); // Mersenne Twister PRNG
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen);
}