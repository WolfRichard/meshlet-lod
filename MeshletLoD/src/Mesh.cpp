#include "Mesh.h"

#include <deque>
#include <metis.h>
#include <unordered_set>

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

    m_meshlets.push_back(dummyMeshlet);

    for (uint i = 0; i < meshlet_count; i++)
    {
        S_Meshlet newMeshlet;
        newMeshlet.triangle_count = meshlets[i].triangle_count;
        newMeshlet.triangle_offset = meshlets[i].triangle_offset;
        newMeshlet.vertex_count = meshlets[i].vertex_count;
        newMeshlet.vertex_offset = meshlets[i].vertex_offset;


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
        S_Meshlet meshlet = m_meshlets[meshlet_index];
        meshlet_edges.push_back(extractEdges(meshlet));
    }
    std::vector<std::vector<uint>> connectivity_matrix;

    for (uint a = 0; a < m_current_hierarchy_top_level_meshlets.size(); a++)
    {
        uint meshlet_A_index = m_current_hierarchy_top_level_meshlets[a];
        S_Meshlet meshlet_A = m_meshlets[meshlet_A_index];
        connectivity_matrix.push_back(std::vector<uint>());
        MATIS_nodeAdjacencyOffsets.push_back((int)MATIS_adjacencyList.size());
        for (uint b = 0; b < m_current_hierarchy_top_level_meshlets.size(); b++)
        {
            uint meshlet_B_index = m_current_hierarchy_top_level_meshlets[b];
            S_Meshlet meshlet_B = m_meshlets[meshlet_B_index];
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
        current_group.hierarchy_tree_depth = m_hierarchy_per_level_group_count.size();

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
            /*
            current_group.bounding_sphere.center.x += current_meshlet.bounding_sphere.center.x;
            current_group.bounding_sphere.center.y += current_meshlet.bounding_sphere.center.y;
            current_group.bounding_sphere.center.z += current_meshlet.bounding_sphere.center.z;
            */
            meshlet_bounding_spheres.push_back(current_meshlet.bounding_sphere);

            for (uint i = 0; i < current_meshlet.triangle_count * 3; i++) 
            {
                merged_deduplicated_indices.push_back(m_vertex_indices[current_meshlet.vertex_offset + m_primitive_indices[current_meshlet.triangle_offset + i]]);
            }
        }
        // average out the bounding sphere centers of chid meshlets to get approximateion for group center
        /*
        current_group.bounding_sphere.center.x /= current_group.meshlet_count;
        current_group.bounding_sphere.center.y /= current_group.meshlet_count;
        current_group.bounding_sphere.center.z /= current_group.meshlet_count;
        */
        current_group.bounding_sphere = computeBoundingSphereRitter(meshlet_bounding_spheres);

        // simplify the meshlet group to have enough space to hold 2 meshlets
        size_t target_index_count = GROUP_SPLIT_COUNT * MAX_MESHLET_VERTEX_COUNT * 3; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should be "GROUP_SPLIT_COUNT * MAX_MESHLET_PRIMITIVE_COUNT * 3", but that would blow past the vertex limit and generate to many meshlets
        std::vector<uint> simplified_indices(merged_deduplicated_indices.size());
        float lod_error = 0.0f;
        size_t simplifiedIndexCount = meshopt_simplify(
            simplified_indices.data(), merged_deduplicated_indices.data(), merged_deduplicated_indices.size(),
            reinterpret_cast<const float*>(&m_vertices.data()[0].position.x), m_vertices.size(), sizeof(S_Vertex),
            target_index_count, FLT_MAX, meshopt_SimplifyLockBorder & meshopt_SimplifyErrorAbsolute, &lod_error);
        simplified_indices.resize(simplifiedIndexCount);
        //current_group.bounding_sphere.radius = lod_error;
    
        

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

            meshopt_Bounds bounds = meshopt_computeMeshletBounds(&m_vertex_indices[newMeshlet.vertex_offset],
                &m_primitive_indices[newMeshlet.triangle_offset],
                newMeshlet.triangle_count,
                &m_vertices[0].position.x,
                m_vertices.size(),
                sizeof(S_Vertex)
            );
            newMeshlet.bounding_sphere.center = float3(bounds.center[0], bounds.center[1], bounds.center[2]);
            newMeshlet.bounding_sphere.radius = bounds.radius;

            m_current_hierarchy_top_level_meshlets.push_back((uint)m_meshlets.size());

            current_group.simplified_meshlets[i] = (uint)m_meshlets.size();

            m_meshlets.push_back(newMeshlet);
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
    current_group.hierarchy_tree_depth = m_hierarchy_per_level_group_count.size();


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

S_BoundingSphere Mesh::computeBoundingSphereRitter(std::vector<S_BoundingSphere>& bounding_spheres) 
{
    S_BoundingSphere result;
    result.center = float3(0, 0, 0);
    result.radius = 0;
    if (bounding_spheres.empty()) return result;

    // Step 1: Find the two most distant sphere centers
    XMVECTOR p1 = XMLoadFloat3(&bounding_spheres[0].center);
    XMVECTOR p2 = XMLoadFloat3(&bounding_spheres[0].center);
    float maxDist = 0;

    for (const auto& s1 : bounding_spheres) {
        for (const auto& s2 : bounding_spheres) {
            float dist = XMVectorGetX(XMVector3LengthSq(XMVectorSubtract(p1, p2)));
            if (dist > maxDist) {
                maxDist = dist;
                p1 = XMLoadFloat3(&s1.center);
                p2 = XMLoadFloat3(&s2.center);
            }
        }
    }

    // Step 2: Compute initial bounding sphere
    XMVECTOR center = (p1 + p2) * 0.5f;
    float radius = XMVectorGetX(XMVector3Length(XMVectorSubtract(p1, p2))) * 0.5f;

    // Step 3: Expand the sphere to include all spheres
    for (const auto& sphere : bounding_spheres) {
        float distToCenter = XMVectorGetX(XMVector3Length(XMVectorSubtract(XMLoadFloat3(&sphere.center), center)));
        if (distToCenter + sphere.radius > radius) {
            float newRadius = (radius + (distToCenter + sphere.radius)) * 0.5f;
            float shift = newRadius - radius;
            center = XMVectorAdd(center, XMVectorScale(XMVector3Normalize(XMVectorSubtract(XMLoadFloat3(&sphere.center), center)), shift));
            radius = newRadius;
        }
    }

    XMStoreFloat3(&result.center, center);
    result.radius = radius;

    return result;
}
