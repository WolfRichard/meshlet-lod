#pragma once

#include <meshoptimizer.h>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include <chrono>
#include <map>
#include <unordered_set>
#include <utility>  // For std::pair
#include <functional>  // For std::hash

struct PairHash {
    // operator() should return a hash value based on the pair
    std::size_t operator()(const std::pair<uint, uint>& p) const {
        // Hash the first and second elements of the pair
        return std::hash<uint>()(p.first) ^ (std::hash<uint>()(p.second) << 1);
    }
};



namespace
{

#include "../shaders/ViewDependentStructures.fxh"

} // namespace

class Mesh
{
public:
    Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);


    std::vector<S_Vertex>        m_vertices;
    std::vector<uint>            m_vertex_indices;
    std::vector<unsigned char>   m_primitive_indices;
    std::vector<uint>            m_morph_indices;
    std::vector<S_Meshlet>       m_meshlets;
    std::vector<S_MeshletGroup>  m_meshlet_groups;
    uint                         m_hierarchy_root_group;
    std::vector<uint>            m_hierarchy_per_level_group_count;

    S_BoundingSphere m_bounding_sphere;
    std::vector<uint> m_original_indices; // indices of the mesh before simplification and meshlet generation


private:
    bool useCustomSimplification = true;


    std::pair<uint, uint> sortEdgeIndices(uint v0, uint v1);
    std::unordered_set<std::pair<uint, uint>, PairHash> extractEdges(S_Meshlet meshlet);
    uint count_shared_edges(const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesA, const std::unordered_set<std::pair<uint, uint>, PairHash>& edgesB);
    std::vector<std::pair<uint, uint>> extractBoundaryEdges(std::vector<uint>& indices);
    
    int randomInt(int min, int max);

    void parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);
    void generateLeafMeshlets();
    void buildMeshletHierachy();
    void groupMeshlets();
    void simplifiyTopLevelGroups();
    void finalTopLevelMeshletGrouping();
    void findParentsItterative();
    S_MeshletGroup getDefaultMeshletGroup();
    void printTreeFromTop(uint current_group_index);
    void printAllMeshlets();
    S_BoundingSphere computeGroupBoundingSphere(S_MeshletGroup& meshlet_group);


    std::vector<uint> m_current_hierarchy_top_level_meshlets; // indices to the meshlets that are still beeing processed in the construction of the DAG
    std::vector<uint> m_current_hierarchy_top_level_groups; // indices to the groups that are still beeing processed in the construction of the DAG
};
