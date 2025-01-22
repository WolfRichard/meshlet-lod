#include "Animation.h"

using namespace DirectX;

void Animation::init(const std::string& animationPath, MeshletMesh* mesh)
{
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(animationPath, aiProcess_Triangulate | aiProcess_MakeLeftHanded | aiProcess_FlipWindingOrder);
    assert(scene && scene->mRootNode);
    auto animation = scene->mAnimations[0];
    m_Duration = (float)animation->mDuration;
    m_TicksPerSecond = (float)animation->mTicksPerSecond;
    ReadHeirarchyData(m_RootNode, scene->mRootNode);
    ReadMissingBones(animation, *mesh);
}

Bone* Animation::FindBone(const std::string& name)
{
    auto iter = std::find_if(m_Bones.begin(), m_Bones.end(),
        [&](const Bone& Bone)
        {
            return Bone.m_Name == name;
        }
    );
    if (iter == m_Bones.end()) return nullptr;
    else return &(*iter);
}

void Animation::ReadMissingBones(const aiAnimation* animation, MeshletMesh& mesh)
{
    int size = animation->mNumChannels;

    auto& boneInfoMap = mesh.m_BoneInfoMap;//getting m_BoneInfoMap from Model class
    int& boneCount = mesh.m_BoneCounter; //getting the m_BoneCounter from Model class

    //reading channels(bones engaged in an animation and their keyframes)
    for (int i = 0; i < size; i++)
    {
        auto channel = animation->mChannels[i];
        std::string boneName = channel->mNodeName.data;

        if (boneInfoMap.find(boneName) == boneInfoMap.end())
        {
            boneInfoMap[boneName].id = boneCount;
            boneCount++;
        }
        m_Bones.push_back(Bone(channel->mNodeName.data,
            boneInfoMap[channel->mNodeName.data].id, channel));
    }

    m_BoneInfoMap = boneInfoMap;
}

void Animation::ReadHeirarchyData(AssimpNodeData& dest, const aiNode* src)
{
    assert(src);

    dest.name = src->mName.data;
    dest.transformation = MeshletMesh::assimpToFloat4x4Matrix(src->mTransformation);
    dest.childrenCount = src->mNumChildren;

    for (uint i = 0; i < src->mNumChildren; i++)
    {
        AssimpNodeData newData;
        ReadHeirarchyData(newData, src->mChildren[i]);
        dest.children.push_back(newData);
    }
}