#pragma once

#include <meshoptimizer.h>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"


namespace
{

#include "../shaders/Structures.fxh"

} // namespace

class MeshletMesh
{
public:
    MeshletMesh(aiMesh* assimp_mesh);

    std::vector<CustomVertex>   m_vertices;
    std::vector<unsigned int>   m_meshlet_vertices;
    std::vector<unsigned char>  m_meshlet_triangles;
    std::vector<DrawTask>       m_draw_tasks;

    float3                      m_boundingSphereCentre;
    float                       m_boundingSphereRadius;

    size_t                      m_index_count;
    size_t                      m_vertex_count;
private:
    void parseMesh(aiMesh* assimp_mesh);
    void generateMeshlets();


    
    std::vector<unsigned int> m_indices;
    

    
    
};
