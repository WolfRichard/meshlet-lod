#pragma once

#include <meshoptimizer.h>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include <chrono>
#include <map>

namespace
{

#include "../shaders/Structures.fxh"

} // namespace

struct BoneInfo
{
    int id; // relative position in its mesh's boneMatrices buffer
    float4x4 transform_matrix; // transformation matrix from model to bone space
};

struct PreBakedAnimation
{
    uint totalAnimationIndex;
    uint boneCount;
    uint frameCount;
    std::vector<float4x4> matrices;
};

class MeshletMesh
{
public:
    MeshletMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);
    static float4x4 assimpToFloat4x4Matrix(const aiMatrix4x4& assimpMatrix)
    {
        return XMMatrixTranspose(DirectX::XMMATRIX(
            assimpMatrix.a1, assimpMatrix.a2, assimpMatrix.a3, assimpMatrix.a4,  
            assimpMatrix.b1, assimpMatrix.b2, assimpMatrix.b3, assimpMatrix.b4,  
            assimpMatrix.c1, assimpMatrix.c2, assimpMatrix.c3, assimpMatrix.c4,  
            assimpMatrix.d1, assimpMatrix.d2, assimpMatrix.d3, assimpMatrix.d4   
        ));
    }

    std::vector<CustomVertex>   m_vertices;
    std::vector<std::vector<unsigned int>>   m_meshlet_vertices;
    std::vector<std::vector<unsigned char>>  m_meshlet_triangles;
    std::vector<std::vector<DrawTask>>       m_draw_tasks;

    float3                      m_boundingSphereCentre;
    float                       m_boundingSphereRadius;

    size_t                      m_index_count;
    size_t                      m_vertex_count;

    std::chrono::duration<double> m_LoDGenTime;
    std::chrono::duration<double> m_meshletGenTime;

    std::map<std::string, BoneInfo> m_BoneInfoMap;
    int m_BoneCounter = 0;

    std::vector<PreBakedAnimation> m_animations;

    std::vector<std::string> m_animationNames;
    std::vector<const char*> m_animationNamesCharP;

    uint m_LoDCount;

private:
    void parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);
    void generateMeshlets();
    void generateLoD(uint resolution_level);

    std::vector<unsigned int> m_indices;
    std::vector<unsigned int> m_simplifiedIndices;

    void addVertexBoneData(CustomVertex& vertex, int boneID, float weight);
    void extractAssimpBoneData(std::vector<CustomVertex>& vertices, aiMesh* mesh, const aiScene* scene);
};
