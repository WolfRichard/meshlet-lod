#include "Mesh.h"

#include <deque>
#include <metis.h>
#include <unordered_set>

#include <random>

#include "Windows.h"

using namespace DirectX;

Mesh::Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    parseMesh(assimp_mesh, assimp_scene);
    generateLeafMeshlets();
    buildMeshletHierachy();

    printTreeFromTop(m_hierarchy_root_group);
    printAllMeshlets();
}


void Mesh::parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
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

        // calc bounding sphere 
        m_bounding_sphere.center.x += new_vert.position.x;
        m_bounding_sphere.center.y += new_vert.position.y;
        m_bounding_sphere.center.z += new_vert.position.z;
    }
    m_bounding_sphere.center.x /= (uint)raw_vertices.size();
    m_bounding_sphere.center.y /= (uint)raw_vertices.size();
    m_bounding_sphere.center.z /= (uint)raw_vertices.size();

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
}

void Mesh::generateLeafMeshlets()
{
    const float cone_weight = 0.0f;
    size_t      max_meshlets = meshopt_buildMeshletsBound(m_original_indices.size(), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT);
    std::vector<meshopt_Meshlet> meshlets(max_meshlets);
    m_vertex_indices.resize(max_meshlets * MAX_MESHLET_VERTEX_COUNT);
    m_primitive_indices.resize(max_meshlets * MAX_MESHLET_PRIMITIVE_COUNT * 3);

    uint meshlet_count = (uint)meshopt_buildMeshlets(meshlets.data(), m_vertex_indices.data(), m_primitive_indices.data(), m_original_indices.data(), m_original_indices.size(), &m_vertices[0].position.x, m_vertices.size(), sizeof(S_Vertex), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT, cone_weight);
    meshlets.resize(meshlet_count);

    S_Meshlet dummyMeshlet; // used to balance holes in groups that cant be filled up to 4 meshlets, always at id 0
    dummyMeshlet.triangle_count = 0;
    dummyMeshlet.triangle_offset = 0;
    dummyMeshlet.vertex_count = 0;
    dummyMeshlet.vertex_offset = 0;
    dummyMeshlet.base_error = 0;
    dummyMeshlet.simplification_error = FLT_MAX;
    dummyMeshlet.bounding_sphere = { float3(0, 0, 0), 0 };

    m_meshlets.push_back(dummyMeshlet);

    for (uint i = 0; i < meshlet_count; i++)
    {
        S_Meshlet newMeshlet;
        newMeshlet.triangle_count = meshlets[i].triangle_count;
        newMeshlet.triangle_offset = meshlets[i].triangle_offset;
        newMeshlet.vertex_count = meshlets[i].vertex_count;
        newMeshlet.vertex_offset = meshlets[i].vertex_offset;
        newMeshlet.base_error = 0;
        newMeshlet.simplification_error = FLT_MAX;
        newMeshlet.simplified_group_bounds = { float3(0, 0, 0), 0 };

        meshopt_Bounds bounds = meshopt_computeMeshletBounds(&m_vertex_indices[newMeshlet.vertex_offset],
            &m_primitive_indices[newMeshlet.triangle_offset],
            newMeshlet.triangle_count,
            &m_vertices[0].position.x,
            m_vertices.size(),
            sizeof(S_Vertex)
        );
        newMeshlet.bounding_sphere.center = float3(bounds.center[0], bounds.center[1], bounds.center[2]);
        newMeshlet.bounding_sphere.radius = bounds.radius;

        m_current_hierarchy_top_level_meshlets.push_back((int)m_meshlets.size());
        m_meshlets.push_back(newMeshlet);
    }
}

void Mesh::buildMeshletHierachy() 
{
    while ((uint)m_current_hierarchy_top_level_meshlets.size() >= GROUP_MERGE_COUNT * 2)
    {
        groupMeshlets();
        simplifiyTopLevelGroups();
    }
    
    assert(m_current_hierarchy_top_level_groups.size() == 2);
    OutputDebugString(("\n\n\nNumber of Top Level Meshlets at tree top: " + std::to_string(m_current_hierarchy_top_level_meshlets.size())).c_str());
    finalTopLevelMeshletGrouping();
    simplifiyTopLevelGroups();

    //clear uneeded data
    m_current_hierarchy_top_level_groups.clear();
    m_current_hierarchy_top_level_meshlets.clear();


    // Check if one of a groups simplfied Meshlets is part of another groups base meshlets, if so declare it its parent
    // with the limitation of only one parent and two children per group

    OutputDebugString(("\n\nNumber of Groups: " + std::to_string(m_meshlet_groups.size()) + "\n\n").c_str());
    OutputDebugString(("\nRecursion Tree Depth: " + std::to_string(m_hierarchy_per_level_group_count.size()) + "\n\n").c_str());
    for (uint count : m_hierarchy_per_level_group_count) OutputDebugString((std::to_string(count) + ", ").c_str());
    findParentsItterative();
    OutputDebugString("\n\n");

    
    for (S_MeshletGroup& current_group : m_meshlet_groups)
    {
        for (uint m = 0; m < current_group.meshlet_count; m++)
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
            current_meshlet.discrete_level_of_detail = current_group.hierarchy_tree_depth;
        }
    }
    

    OutputDebugString(("Size vertex_indices vector: " + std::to_string(m_vertex_indices.size()) + " Size morph_indices vector: " + std::to_string(m_morph_indices.size()) + "\n").c_str());
}



void Mesh::groupMeshlets()
{
    //while ((int)m_current_hierarchy_top_level_meshlets.size() % GROUP_MERGE_COUNT != 0)
    //    m_current_hierarchy_top_level_meshlets.push_back(0); // fill with dummy meshlets until equal meshlet distribution is possible between all groups

    // build METIS graph representation of meshlet connectivity to each other
    idx_t MATIS_numNodes = (int)m_current_hierarchy_top_level_meshlets.size();
    std::vector<idx_t> MATIS_vertexWeights(MATIS_numNodes, 100);
    idx_t MATIS_numEdges = 0;
    std::vector<idx_t> MATIS_nodeAdjacencyOffsets;
    std::vector<idx_t> MATIS_adjacencyList;
    std::vector<idx_t> MATIS_edgeWeights; // represents the number of shared edges between each meshlet

    // extract all edges from the current top level meshlets
    std::vector<std::unordered_set<std::pair<uint, uint>, PairHash>> meshlet_edges;
    for (auto meshlet_index : m_current_hierarchy_top_level_meshlets)
    {
        meshlet_edges.push_back(extractEdges(m_meshlets[meshlet_index]));
    }
    std::vector<std::vector<uint>> connectivity_matrix;

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

    
    OutputDebugString("\n\n\n");

    OutputDebugString(("Number of Meshlets: " + std::to_string(m_current_hierarchy_top_level_meshlets.size())).c_str());

    int debugCounter = 0;
    for (auto& meshletSharedEdgesList : connectivity_matrix)
    {
        OutputDebugString(("\nMeshlet_" + std::to_string(debugCounter++) + ": ").c_str());
        uint neighboringMeshletCount = 0;
        for (auto matrix_entry : meshletSharedEdgesList)
        {
            OutputDebugString((std::to_string(matrix_entry) + ",").c_str());
            if (matrix_entry) neighboringMeshletCount++;
        }
        if (neighboringMeshletCount) OutputDebugString(("\nMeshlet has " + std::to_string(neighboringMeshletCount) + " Neighbors").c_str());
        else OutputDebugString("\nMeshlet shares no Edge!!!!!");

    }

    OutputDebugString("\nRESULT: ");
    for (auto assignedGroupIndex : MATIS_outputPartitionsArray)
        OutputDebugString((std::to_string(assignedGroupIndex) + ", ").c_str());
    OutputDebugString("\nNumber of Meshlets in each group: ");
    for (int k = 0; k < MATIS_numPartitions; k++) 
    {
        debugCounter = 0;
        for (auto assignedGroupIndex : MATIS_outputPartitionsArray)
            if (k == assignedGroupIndex)
                debugCounter++;
        OutputDebugString((std::to_string(debugCounter) + ", ").c_str());
    }

    OutputDebugString("\n\n\n");
    

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


// Sorts indices of an edge from small to big
std::pair<uint, uint> Mesh::sortEdgeIndices(uint v0, uint v1) 
{
    return (v0 < v1) ? std::make_pair(v0, v1) : std::make_pair(v1, v0);
}


// Function to extract edges from a meshlet
std::unordered_set<std::pair<uint, uint>, PairHash> Mesh::extractEdges(S_Meshlet meshlet) 
{
    std::unordered_set<std::pair<uint, uint>, PairHash> edge_set;

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

// Function to count shared edges between two meshlets
uint Mesh::count_shared_edges(const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesA, const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesB) 
{
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

    // Step 1: Count occurrences of each edge
    for (uint tri = 0; tri < (uint)indices.size() / 3; tri++)
    {
        edgeCount[sortEdgeIndices(indices[tri * 3 + 0], indices[tri * 3 + 1])]++;
        edgeCount[sortEdgeIndices(indices[tri * 3 + 1], indices[tri * 3 + 2])]++;
        edgeCount[sortEdgeIndices(indices[tri * 3 + 2], indices[tri * 3 + 0])]++;
    }

    // Step 2: Collect edges that appear only once (boundary edges)
    std::vector<std::pair<uint, uint>> boundaryEdges;
    for (const auto& [edge, count] : edgeCount) {
        if (count == 1) {  // Boundary edges appear only once
            boundaryEdges.push_back(edge);
        }
    }

    return boundaryEdges;
}

void Mesh::simplifiyTopLevelGroups()
{
    m_current_hierarchy_top_level_meshlets.clear(); // clear the top level of meshlets as it will be replaced by their simplified version

    OutputDebugString(("\n\nNumber of Top LEVEL GROUPS: " + std::to_string(m_current_hierarchy_top_level_groups.size())).c_str());
    for (uint g = 0; g < m_current_hierarchy_top_level_groups.size(); g++) 
    {
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
         
        computeGroupBoundingSphere(current_group);
        //current_group.bounding_sphere = computeBoundingSphereRitter(meshlet_bounding_spheres);

        // simplify the meshlet group to have enough space to hold 2 meshlets
        size_t target_index_count = GROUP_SPLIT_COUNT * MAX_MESHLET_VERTEX_COUNT * 3; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should be "GROUP_SPLIT_COUNT * MAX_MESHLET_PRIMITIVE_COUNT * 3", but that would blow past the vertex limit and generate to many meshlets
        std::vector<uint> simplified_indices(merged_deduplicated_indices.size());
        float lod_error = 0.0f;
        size_t simplifiedIndexCount = meshopt_simplify(
            simplified_indices.data(), merged_deduplicated_indices.data(), merged_deduplicated_indices.size(),
            reinterpret_cast<const float*>(&m_vertices.data()[0].position.x), m_vertices.size(), sizeof(S_Vertex),
            target_index_count, FLT_MAX, meshopt_SimplifyLockBorder, &lod_error); // TODO: & meshopt_SimplifyErrorAbsolute bitmask option leads to cracks??? why??
        simplified_indices.resize(simplifiedIndexCount);
        //current_group.bounding_sphere.radius = lod_error;













        // Custom simplification algorithm that tracks vertex decimation indices

        if (useCustomSimplification) 
        {   
            struct simpl_triangle
            {
                uint vertex_indices[3]; // index into the simplification vertex struct array
                uint edge_indices[3]; // index into the simplification edge struct array
                bool removed = false; // true if element was removed during simplification process
            };

            struct simpl_edge
            {
                std::pair<uint, uint> vertex_indices; // index into the simplification vertex struct array
                std::vector<uint> triangle_indices; // edge can be part of 1 or 2 triangles
                bool is_boundary_edge; // true if edge forms the boundary of the geometry
                bool removed = false; // true if element was removed during simplification process
            };

            struct simpl_vertex
            {
                std::vector<uint> edge_indices; // edges that connect to the vertex
                std::vector<uint> triangle_indices; // triangles that connect to the vertex
                uint vertex_index; // index into the m_vertices vector
                std::vector<uint> merged_vertex_indices; // indices that point into the s_vertex struct array that have been merged into this vertex
                bool is_boundary_vertex = false; // true if vertex lies on the edge of the geometry
                bool removed = false; // true if element was removed during simplification process
            };

            std::vector<simpl_vertex> simplification_vertices;
            std::vector<simpl_edge> simplification_edges;
            std::vector<simpl_triangle> simplification_triangles;

          
            // go over every triangle of the original index buffer
            for (uint t = 0; t < merged_deduplicated_indices.size() / 3; t++)
            {
                // check if ths is a degenerated triangle, if so skip it
                if (   merged_deduplicated_indices[t * 3 + 0] == merged_deduplicated_indices[t * 3 + 1]
                    || merged_deduplicated_indices[t * 3 + 1] == merged_deduplicated_indices[t * 3 + 2]
                    || merged_deduplicated_indices[t * 3 + 2] == merged_deduplicated_indices[t * 3 + 0])
                {
                    OutputDebugString("DEGENERATE TRIANGLE!!!!\n");
                    assert(false);
                    continue;
                }

                simpl_triangle new_triangle;
                // search for each of the vertices structs that are part of this triangles. create if vertex was not yet created
                for (uint i = 0; i < 3; i++)
                {
                    int s_vertex_index = -1;
                    for (uint svi = 0; svi < simplification_vertices.size(); svi++)
                    {
                        simpl_vertex& sv = simplification_vertices[svi];
                        if (sv.vertex_index == merged_deduplicated_indices[t * 3 + i])
                        {
                            s_vertex_index = svi;
                            break;
                        }
                    }
                    if (s_vertex_index == -1)
                    {
                        s_vertex_index = simplification_vertices.size();
                        simpl_vertex new_vertex;
                        new_vertex.vertex_index = merged_deduplicated_indices[t * 3 + i];
                        simplification_vertices.push_back(new_vertex);
                    }
                    new_triangle.vertex_indices[i] = s_vertex_index;
                    simplification_vertices[s_vertex_index].triangle_indices.push_back(simplification_triangles.size());
                }

                // search and if not already present create a new edge struct for each of the triangle sides
                for (uint i = 0; i < 3; i++) 
                {
                    simpl_edge new_edge;
                    new_edge.vertex_indices = sortEdgeIndices(new_triangle.vertex_indices[i], new_triangle.vertex_indices[(i + 1) % 3]);
                    int s_edge_index = -1;
                    for (uint sei = 0; sei < simplification_edges.size(); sei++)
                    {
                        simpl_edge& se = simplification_edges[sei];
                        if (se.vertex_indices == new_edge.vertex_indices)
                        {
                            s_edge_index = sei;
                            break;
                        }
                    }
                    if (s_edge_index == -1)
                    {
                        s_edge_index = simplification_edges.size();
                        simplification_edges.push_back(new_edge);
                    }
                    new_triangle.edge_indices[i] = s_edge_index;
                    assert(new_edge.vertex_indices.first != new_edge.vertex_indices.second);
                    simplification_edges[s_edge_index].triangle_indices.push_back(simplification_triangles.size());
                }
                simplification_triangles.push_back(new_triangle);
            }
               
            // go over every edge, update boundary information for each edge and vertex while also adding the connectivity data from vertex to edge
            for (uint sei = 0; sei < simplification_edges.size(); sei++)
            {
                simpl_edge& se = simplification_edges[sei];
                //OutputDebugString(("Number of adjecent Triangles to edge: " + std::to_string(se.triangle_indices.size()) + "\n").c_str());
                if (se.triangle_indices.size() > 2)
                {
                    OutputDebugString(("Following Edge has to many triangles: " + std::to_string(se.vertex_indices.first) + ", " + std::to_string(se.vertex_indices.second) + "\n").c_str());
                    for (uint t : se.triangle_indices)
                    {
                        simpl_triangle& st = simplification_triangles[t];
                        OutputDebugString(("Triangles adjecent to this edge: " + std::to_string(st.vertex_indices[0]) + ", " + std::to_string(st.vertex_indices[1]) + ", " + std::to_string(st.vertex_indices[2]) + "\n").c_str());
                    }
                }
                assert((se.triangle_indices.size() >= 1) && (se.triangle_indices.size() <= 2)); // TODO: can have more than 2 triangles, because of differing winding order so technically incorrect!
                se.is_boundary_edge = (se.triangle_indices.size() == 1); // edge forms boundary if its only part of one triangle

                // add connectivity data from vertex to edge
                simplification_vertices[se.vertex_indices.first].edge_indices.push_back(sei);
                simplification_vertices[se.vertex_indices.second].edge_indices.push_back(sei);

                // inherit boundary satus from edges to their vertices
                simplification_vertices[se.vertex_indices.first].is_boundary_vertex = se.is_boundary_edge;
                simplification_vertices[se.vertex_indices.second].is_boundary_vertex = se.is_boundary_edge;
            }

            // extract vertex indices of all non bondary simpl_vertices
            std::vector<uint> vertex_valid_for_simplification_pool;
            for (uint i = 0; i < simplification_vertices.size(); i++)
            {
                if (!(simplification_vertices[i].is_boundary_vertex))
                    vertex_valid_for_simplification_pool.push_back(i);
            }


            uint current_primitives_count = (uint)simplification_triangles.size();
            uint current_vertex_count = (uint)simplification_vertices.size();

            OutputDebugString(("Triangle Count: " + std::to_string(current_primitives_count) + "\n").c_str());
            OutputDebugString(("Vertex Count: " + std::to_string(current_vertex_count) + "\n").c_str());
            OutputDebugString(("Pool Vertex Count: " + std::to_string(vertex_valid_for_simplification_pool.size()) + "\n").c_str());
            OutputDebugString(("Boundary Vertex Count: " + std::to_string(simplification_vertices.size() - vertex_valid_for_simplification_pool.size()) + "\n").c_str());

            // itterativly decimate vertices
            while ((current_primitives_count > MAX_MESHLET_PRIMITIVE_COUNT * GROUP_SPLIT_COUNT) 
                || (current_vertex_count > MAX_MESHLET_VERTEX_COUNT))// * GROUP_SPLIT_COUNT)) // TODO: find a better solu´tion to the problem
            {
                // pick vertex that should be removed
                // TODO: currently chooses random vertex from all inner vertices that have the smallest number of edges connected to them. this should be replaced by QEM in the future!
                assert(vertex_valid_for_simplification_pool.size() > 0);

                std::vector<uint> lowest_connectivity_vertex_index_pool_indices;

                uint lowest_edge_count = UINT32_MAX;
                for (uint vvfsi = 0; vvfsi < vertex_valid_for_simplification_pool.size(); vvfsi++)
                {
                    simpl_vertex& sv = simplification_vertices[vertex_valid_for_simplification_pool[vvfsi]];
                    
                    if (sv.edge_indices.size() == lowest_edge_count)
                    {
                        lowest_connectivity_vertex_index_pool_indices.push_back(vvfsi);
                    }
                    else if (sv.edge_indices.size() < lowest_edge_count)
                    {
                        lowest_edge_count = (uint)sv.edge_indices.size();
                        lowest_connectivity_vertex_index_pool_indices.clear();
                        lowest_connectivity_vertex_index_pool_indices.push_back(vvfsi);
                    }
                    
                }


                uint randomly_selected_vertex_simplification_pool_index = lowest_connectivity_vertex_index_pool_indices[randomInt(0, lowest_connectivity_vertex_index_pool_indices.size() - 1)];
                uint removing_vertex_index = vertex_valid_for_simplification_pool[randomly_selected_vertex_simplification_pool_index];
                simpl_vertex& removing_vertex = simplification_vertices[removing_vertex_index];
                S_Vertex& dereferenced_removing_vertex = m_vertices[removing_vertex.vertex_index];

                // find closest neighbor vertex that should be selected as the morph target
                float min_dist2 = FLT_MAX;
                uint closest_vertex_neighbor_index;
                uint shortest_neighbor_edge_index;
                for (uint ei : removing_vertex.edge_indices)
                {
                    simpl_edge& e = simplification_edges[ei];
                    assert(e.vertex_indices.first != e.vertex_indices.second);
                    uint other_edge_vertex_index = (e.vertex_indices.first == removing_vertex_index) ? e.vertex_indices.second : e.vertex_indices.first;
                    simpl_vertex& other_vertex = simplification_vertices[other_edge_vertex_index];
                    S_Vertex dereferenced_other_vertex = m_vertices[other_vertex.vertex_index];

                    float dx = dereferenced_removing_vertex.position.x - dereferenced_other_vertex.position.x;
                    float dy = dereferenced_removing_vertex.position.y - dereferenced_other_vertex.position.y;
                    float dz = dereferenced_removing_vertex.position.z - dereferenced_other_vertex.position.z;
                    float dist2 = dx * dx + dy * dy + dz * dz;

                    if (dist2 < min_dist2)
                    {
                        shortest_neighbor_edge_index = ei;
                        closest_vertex_neighbor_index = other_edge_vertex_index;
                        min_dist2 = dist2;
                    }
                }


                
                // remove flag the selected vertex, its shortest neighboring edge and all triangles that side the removed edge
                simpl_edge& removing_edge = simplification_edges[shortest_neighbor_edge_index];
                for (uint removing_triangle_index : removing_edge.triangle_indices)
                {
                    simpl_triangle& removing_triangle = simplification_triangles[removing_triangle_index];
                    removing_triangle.removed = true;
                    current_primitives_count--;
                }
                removing_edge.removed = true;
                removing_vertex.removed = true;
                assert((sortEdgeIndices(closest_vertex_neighbor_index, removing_vertex_index) == removing_edge.vertex_indices));
                current_vertex_count--;

                
                assert(closest_vertex_neighbor_index != removing_vertex_index);
                // append all morph children vertex indices into own morph target aswell as removed vertex index
                simpl_vertex& morph_target_vertex = simplification_vertices[closest_vertex_neighbor_index];
                morph_target_vertex.merged_vertex_indices.insert(morph_target_vertex.merged_vertex_indices.end(), removing_vertex.merged_vertex_indices.begin(), removing_vertex.merged_vertex_indices.end());
                morph_target_vertex.merged_vertex_indices.push_back(removing_vertex_index);

                // erase old connection indices from removed edge and erased triangles from morph target
                for (int morph_target_triangles_index = morph_target_vertex.triangle_indices.size() - 1; morph_target_triangles_index >= 0; morph_target_triangles_index--)
                {
                    if (simplification_triangles[morph_target_vertex.triangle_indices[morph_target_triangles_index]].removed)
                    {
                        morph_target_vertex.triangle_indices[morph_target_triangles_index] = morph_target_vertex.triangle_indices.back();
                        morph_target_vertex.triangle_indices.pop_back();
                    }
                }
                for (int morph_target_edge_index = morph_target_vertex.edge_indices.size() - 1; morph_target_edge_index >= 0; morph_target_edge_index--)
                {
                    if (simplification_edges[morph_target_vertex.edge_indices[morph_target_edge_index]].removed)
                    {
                        morph_target_vertex.edge_indices[morph_target_edge_index] = morph_target_vertex.edge_indices.back();
                        morph_target_vertex.edge_indices.pop_back();
                    }
                }

                // update all other neighboring edges and triangles to have their vertex index replaced by the morph target
                // add all other adjacent triangles and edges to the morph target if they wreent removed flagged
                for (uint neighboring_triangle_index : removing_vertex.triangle_indices)
                {
                    simpl_triangle& nt = simplification_triangles[neighboring_triangle_index];
                    if (nt.removed) continue; // skip removed triangles
                    for (uint& vi : nt.vertex_indices)
                    {
                        if (vi == removing_vertex_index) vi = closest_vertex_neighbor_index;
                    }
                    morph_target_vertex.triangle_indices.push_back(neighboring_triangle_index);
                }
                for (uint neighboring_edge_index : removing_vertex.edge_indices)
                {
                    simpl_edge& ne = simplification_edges[neighboring_edge_index];
                    if (ne.removed) continue; // skip removed edge
                    assert(ne.vertex_indices.first != ne.vertex_indices.second);
                    assert(shortest_neighbor_edge_index != neighboring_edge_index);
                    assert(!(sortEdgeIndices(closest_vertex_neighbor_index, removing_vertex_index) == ne.vertex_indices));
                    if (ne.vertex_indices.first == removing_vertex_index)
                    {
                        ne.vertex_indices.first = closest_vertex_neighbor_index;
                    }
                    else
                    {
                        ne.vertex_indices.second = closest_vertex_neighbor_index;
                    }
                    morph_target_vertex.edge_indices.push_back(neighboring_edge_index);
                    assert(ne.vertex_indices.first != ne.vertex_indices.second);
                }

                

                // remove vertex from simplification pool
                vertex_valid_for_simplification_pool[randomly_selected_vertex_simplification_pool_index] = vertex_valid_for_simplification_pool.back();
                vertex_valid_for_simplification_pool.pop_back();
            }
            

            // reconstruct simplified index buffer from custom "simpl_" struct vectors
            simplified_indices.clear();

            for (simpl_triangle& st : simplification_triangles)
            {
                if (!(st.removed))
                {
                    for (uint i = 0; i < 3; i++) {
                        simpl_vertex& sv = simplification_vertices[st.vertex_indices[i]];
                        assert(!(sv.removed));
                        simplified_indices.push_back(sv.vertex_index);
                    }
                }
            }


            OutputDebugString(("Triangle Count after custom simp: " + std::to_string(current_primitives_count) + "\n").c_str());
            OutputDebugString(("Vertex Count after custom simp: " + std::to_string(current_vertex_count) + "\n").c_str());
        }















        
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

        // generate new meshlets on simplified geometry
        const float cone_weight = 0.0f;
        size_t      max_meshlets = meshopt_buildMeshletsBound(simplified_indices.size(), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT);
        std::vector<meshopt_Meshlet> meshlets(max_meshlets);
        std::vector<uint> vertex_indices(max_meshlets * MAX_MESHLET_VERTEX_COUNT);
        std::vector<unsigned char> primitive_indices(max_meshlets * MAX_MESHLET_PRIMITIVE_COUNT * 3);

        uint meshlet_count = (uint)meshopt_buildMeshlets(meshlets.data(), vertex_indices.data(), primitive_indices.data(), simplified_indices.data(), simplified_indices.size(), &m_vertices[0].position.x, m_vertices.size(), sizeof(S_Vertex), MAX_MESHLET_VERTEX_COUNT, MAX_MESHLET_PRIMITIVE_COUNT, cone_weight);
        meshlets.resize(meshlet_count);

        OutputDebugString(("\n\nNumber of Triangles before Simplification: " + std::to_string(merged_deduplicated_indices.size() / 3) + " Number of Triangles after Simplification: " + std::to_string(simplified_indices.size() / 3)).c_str());
        OutputDebugString(("\nNumber of Meshlets: " + std::to_string(meshlet_count) + " Maximum Number of Meshlets: " + std::to_string(max_meshlets)).c_str());
        for (auto meshlet : meshlets)
            OutputDebugString(("\n\tNumber of Vertices: " + std::to_string(meshlet.vertex_count) + " Number of Triangles: " + std::to_string(meshlet.triangle_count)).c_str());

        //assert(meshlet_count == 2); //currently root group only simplifies to a single meshlet???

        // add new vertex/primitive-indicesbuffers to the main buffer and offset meshlet meta data accordingly
        uint previous_primitive_indices_size = (uint)m_primitive_indices.size();
        uint previous_vertex_indices_size = (uint)m_vertex_indices.size();

        m_vertex_indices.insert(m_vertex_indices.end(), vertex_indices.begin(), vertex_indices.end());
        m_primitive_indices.insert(m_primitive_indices.end(), primitive_indices.begin(), primitive_indices.end());

        for (uint i = 0; i < meshlet_count; i++)
        {
            S_Meshlet newMeshlet;
            newMeshlet.triangle_count = meshlets[i].triangle_count;
            newMeshlet.triangle_offset = meshlets[i].triangle_offset + previous_primitive_indices_size;
            newMeshlet.vertex_count = meshlets[i].vertex_count;
            newMeshlet.vertex_offset = meshlets[i].vertex_offset + previous_vertex_indices_size;
            newMeshlet.base_error = collective_base_meshlet_error;
            newMeshlet.simplification_error = FLT_MAX;
            newMeshlet.simplified_group_bounds = { float3(0, 0, 0), 0 };

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

            newMeshlet.bounding_sphere = current_group.bounding_sphere; // !!!!!!!!!!!!!!!!!?

            m_current_hierarchy_top_level_meshlets.push_back((uint)m_meshlets.size());

            current_group.simplified_meshlets[i] = (uint)m_meshlets.size();


            m_meshlets.push_back(newMeshlet);
        }
        
        
        // TODO!!!
        // populate morph_indices buffer data
        //std::vector<std::pair<uint, uint>>  simplified_group_boundary_edges = extractBoundaryEdges(current_group);
        m_morph_indices.resize(m_vertex_indices.size(), 0); // resize in case new meshlets got added in previous hierarchy expansion
        for (uint m = 0; m < current_group.meshlet_count; m++) // go over every meshlet of the current group
        {
            S_Meshlet& current_meshlet = m_meshlets[current_group.meshlets[m]];
         


            // find closest vertex that survived simplification to set as morph target //TODO also morph to boundary vertices from inner ones
            for (uint vi = 0; vi < current_meshlet.vertex_count; vi++) 
            {
                uint vertex_index = m_vertex_indices[current_meshlet.vertex_offset + vi];
                S_Vertex current_vertex = m_vertices[vertex_index];
                float min_dist2 = FLT_MAX;
                uint closest_morph_target_vertex_index_index = current_meshlet.vertex_offset + vi;


                for (uint m = 0; m < GROUP_SPLIT_COUNT; m++) // go over every simplified meshlet of the current group
                {
                    S_Meshlet& current_simplified_meshlet = m_meshlets[current_group.simplified_meshlets[m]];
                    for (unsigned char svi = 0; svi < current_simplified_meshlet.vertex_count; svi++) 
                    {
                        uint simplified_vertex_index = m_vertex_indices[current_simplified_meshlet.vertex_offset + svi];
                        if (vertex_index == simplified_vertex_index)
                        {
                            closest_morph_target_vertex_index_index = current_simplified_meshlet.vertex_offset + svi;
                            min_dist2 = 0;
                            break;
                        }

                        S_Vertex comparing_vertex = m_vertices[simplified_vertex_index];

                        float dx = current_vertex.position.x - comparing_vertex.position.x;
                        float dy = current_vertex.position.y - comparing_vertex.position.y;
                        float dz = current_vertex.position.z - comparing_vertex.position.z;
                        float dist2 = dx * dx + dy * dy + dz * dz;

                        if (dist2 < min_dist2)
                        {
                            min_dist2 = dist2;
                            closest_morph_target_vertex_index_index = current_simplified_meshlet.vertex_offset + svi;
                        }
                    }
                }
                m_morph_indices[current_meshlet.vertex_offset + vi] = closest_morph_target_vertex_index_index;
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


void Mesh::findParentsItterative()
{
    uint accumulated_offset = 0;
    for (uint tree_level = 0; tree_level < m_hierarchy_per_level_group_count.size(); tree_level++)
    {
        OutputDebugString(("Current Tree Level - " + std::to_string(tree_level) + "\n").c_str());

        for (uint i = 0; i < m_hierarchy_per_level_group_count[tree_level]; i++)
        {
            uint current_group_index = accumulated_offset + i;
            S_MeshletGroup& current_group = m_meshlet_groups[current_group_index];
            
            uint possible_parent_with_least_children = -1;
            // search every group to find where the base meshlets contain one of the simplified meshlets of the current group and make it its your parent
            for (uint comparing_group_index = 0; comparing_group_index < m_meshlet_groups.size(); comparing_group_index++)
            {
                S_MeshletGroup& comparing_group = m_meshlet_groups[comparing_group_index];
                
                for (uint comparing_base_meshlet_index = 0; comparing_base_meshlet_index < comparing_group.meshlet_count; comparing_base_meshlet_index++)
                {
                    uint comparing_base_meshlet = comparing_group.meshlets[comparing_base_meshlet_index];
                    for (uint current_simplified_meshlet : current_group.simplified_meshlets)
                    {
                        if (comparing_base_meshlet == current_simplified_meshlet)
                        {
                            if (possible_parent_with_least_children == -1) possible_parent_with_least_children = comparing_group_index;
                            else if (comparing_group.childCount < m_meshlet_groups[possible_parent_with_least_children].childCount) possible_parent_with_least_children = comparing_group_index;
                            
                            //comparing_group.bounding_sphere.radius += current_group.bounding_sphere.radius; // add child error ontop of all parents error to ensure that parents error always stays highter than that of child
                            break; // prevents chlid error to be added multiple times to same parrent if multiple meshlets match with parent
                        }
                    }
                }
            }

            if (possible_parent_with_least_children != -1)
            {
                current_group.parent = possible_parent_with_least_children;
                S_MeshletGroup& parent_group = m_meshlet_groups[possible_parent_with_least_children];
                parent_group.children[parent_group.childCount++] = current_group_index;
                
            }

            OutputDebugString(("Group[" + std::to_string(current_group_index) + "] --> Parent[" + std::to_string(current_group.parent) + "]\n").c_str());
        }
        accumulated_offset += m_hierarchy_per_level_group_count[tree_level];
    }
}

S_MeshletGroup Mesh::getDefaultMeshletGroup()
{
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
    S_MeshletGroup& current_group = m_meshlet_groups[current_group_index];
    OutputDebugString(("\nGroup[" + std::to_string(current_group_index) + "] --> Children: ").c_str());
    for (uint i = 0; i < current_group.childCount; i++)
        OutputDebugString(("[" + std::to_string(current_group.children[i]) + "], ").c_str());
    for (uint i = 0; i < current_group.childCount; i++)
        printTreeFromTop(current_group.children[i]);
}

void Mesh::printAllMeshlets()
{
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

    // estimate center
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