#pragma once

#include "HLSLnames.h"

#include "Bone.h"
#include "MeshletMesh.h"

#include <map>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

struct AssimpNodeData
{
    float4x4 transformation;
    std::string name;
    int childrenCount;
    std::vector<AssimpNodeData> children;
};

class Animation
{
public:
    Animation() {}
    void init(const aiScene* scene, uint animationIndex, MeshletMesh* mesh);
    ~Animation() {}
    Bone* FindBone(const std::string& name);

    float m_Duration = 0;
    float m_TicksPerSecond = 0;
    AssimpNodeData m_RootNode;
    std::map<std::string, BoneInfo> m_BoneInfoMap;
    std::vector<Bone> m_Bones;

private:
    void ReadMissingBones(const aiAnimation* animation, MeshletMesh& mesh);
    void ReadHeirarchyData(AssimpNodeData& dest, const aiNode* src);

    
};