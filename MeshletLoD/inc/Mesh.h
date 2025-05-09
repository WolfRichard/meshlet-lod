#pragma once

#include <meshoptimizer.h>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include "../shaders/ViewDependentStructures.fxh"

#include <chrono>
#include <map>
#include <unordered_set>
#include <utility>  
#include <functional>


// helper function that generates a hash for a given uint pair
// needed for setting up and using a map of pairs when simplifying the mesh
struct PairHash {
    // operator() should return a hash value based on the pair
    std::size_t operator()(const std::pair<uint, uint>& p) const {
        // Hash the first and second elements of the pair
        return std::hash<uint>()(p.first) ^ (std::hash<uint>()(p.second) << 1);
    }
};


// Loads and manages a single unique meshe and its corrisponding data
// automatically generates all needed LOD help structures for later rendering
class Mesh
{
public:
    Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);

    //geometry data
    std::vector<S_Vertex>        m_vertices;
    std::vector<uint>            m_vertex_indices;
    std::vector<unsigned char>   m_primitive_indices;
    std::vector<uint>            m_morph_indices;
    std::vector<S_Meshlet>       m_meshlets;

    // helping data structures for hierarchy group building and simplification
    std::vector<S_MeshletGroup>  m_meshlet_groups;
    uint                         m_hierarchy_root_group;
    std::vector<uint>            m_hierarchy_per_level_group_count;

    // conservative bounding sphere around all vertices of the whole mesh
    S_BoundingSphere m_bounding_sphere;

    // indices of the mesh before simplification and meshlet generation
    std::vector<uint> m_original_indices; 

private:
    // custom simplification needed for geomorphing
    // if disabled it will default to the simplification provided by the mesh optimizer library which is currently faster and results in better simplification quality
    bool m_useCustomSimplification = false;

    // limits the maximum triangle count of meshes before they are preprocessed, so that preprocessing does not exceed reasonable times
    uint m_maximum_mesh_triangle_count = 100000; // 1 million primitives taes around 7 min

    // returns an edge represented as a pair of two uints
    // vertex indices will be ordered by value to make later comparing of two edges simpler
    std::pair<uint, uint> sortEdgeIndices(uint v0, uint v1);
    // extracts all edges of a specified meshlet and stores them inside a hash map for faster look up
    std::unordered_set<std::pair<uint, uint>, PairHash> extractEdges(S_Meshlet meshlet);
    // extracts all boundary edges of a specified meshlet and stores them inside a hash map for faster look up
    std::unordered_set<std::pair<uint, uint>, PairHash> Mesh::extractBoundaryEdgeSet(S_Meshlet meshlet);
    // ccompares the number of boundary edges that are shared between two meshlets
    uint count_shared_edges(const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesA, const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesB);
    // extract all edges that form the border of a meshlet (all edges that are only part of one triangle)
    std::vector<std::pair<uint, uint>> extractBoundaryEdges(std::vector<uint>& indices);
    
    // returns a random integer between and including the specified minimum and maximum values
    int randomInt(int min, int max);

    // loads and parses the data of a single individual unique mesh
    void parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);

    // first step of the hierarchy building
    // generates the base meshlets that represent the original unsimplified mesh geometry
    void generateLeafMeshlets();

    // generates the hierarchical LOD data from the leaf meshlets
    void buildMeshletHierachy();

    // second step in the hierarchy building process, that groups 4 to 5 meshlets by minimizing edge cost through graph partitioning
    // will be called for every level of the hierarchy tree
    void groupMeshlets();

    // simplifies the current level of newly grouped meshlets
    // this function keeps the original meshlets that were simplified and stores the simplification result in new meshlets that will form the foundation for the next tree level
    // each group of meshlets is ittertivly simplified until the new vertex and primitive count is reduced to the size of 2 new meshlets 
    void simplifiyTopLevelGroups();

    // special edge case handeling for the tree root
    void finalTopLevelMeshletGrouping();


    // generates a meshlet group with default initialized values (always conservative values that would result in 0 workload)
    S_MeshletGroup getDefaultMeshletGroup();

    // debugging function
    void printTreeFromTop(uint current_group_index);
    void printAllMeshlets();
    void findParentsItterative();

    // calculates a conservative approximation of the mesh's bounding sphere
    S_BoundingSphere computeGroupBoundingSphere(S_MeshletGroup& meshlet_group);

    // indices to the meshlets that are still beeing processed in the construction of the DAG
    std::vector<uint> m_current_hierarchy_top_level_meshlets; 
    // indices to the groups that are still beeing processed in the construction of the DAG
    std::vector<uint> m_current_hierarchy_top_level_groups; 
};
