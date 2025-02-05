#include "Animator.h"
#include "Windows.h"


using namespace DirectX;

void Animator::recalibrateOffsets(const AssimpNodeData* node, float4x4 parentTransform)
{
    std::string nodeName = node->name;
    float4x4 nodeTransform = node->transformation;

    Bone* Bone = m_CurrentAnimation->FindBone(nodeName);

    if (Bone)
    {
        Bone->Update(0);//m_CurrentTime);
        nodeTransform = Bone->m_LocalTransform;
        //    OutputDebugString((std::string("Did find:") + nodeName + std::string("\n")).c_str());
    }

    float4x4 globalTransformation = nodeTransform * parentTransform;

    auto boneInfoMap = m_CurrentAnimation->m_BoneInfoMap;
    if (boneInfoMap.find(nodeName) != boneInfoMap.end())
    {
        int index = boneInfoMap[nodeName].id;
        float4x4 offset = boneInfoMap[nodeName].transform_matrix;
        m_CurrentAnimation->m_BoneInfoMap[nodeName].transform_matrix = XMMatrixInverse(nullptr, globalTransformation);
    }

    for (int i = 0; i < node->childrenCount; i++)
        recalibrateOffsets(&node->children[i], globalTransformation);
}

void Animator::init(Animation* animation)
{
    m_CurrentTime = 0.0;
    m_CurrentAnimation = animation;
   
    m_FinalBoneMatrices.reserve(animation->m_Bones.size());
    for (int i = 0; i < animation->m_Bones.size(); i++)
        m_FinalBoneMatrices.push_back(DirectX::XMMatrixIdentity());

    //recalibrateOffsets(&m_CurrentAnimation->m_RootNode, DirectX::XMMatrixIdentity());
}

void Animator::UpdateAnimation(float dt)
{
    m_DeltaTime = dt;
    if (m_CurrentAnimation)
    {
        m_CurrentTime += m_CurrentAnimation->m_TicksPerSecond * dt;
        m_CurrentTime = fmod(m_CurrentTime, m_CurrentAnimation->m_Duration);

        CalculateBoneTransform(&m_CurrentAnimation->m_RootNode, DirectX::XMMatrixIdentity());
    }
}

void Animator::PlayAnimation(Animation* pAnimation)
{
    m_CurrentAnimation = pAnimation;
    m_CurrentTime = 0.0f;
}

void Animator::CalculateBoneTransform(const AssimpNodeData* node, float4x4 parentTransform)
{
    std::string nodeName = node->name;
    float4x4 nodeTransform = node->transformation;

    Bone* Bone = m_CurrentAnimation->FindBone(nodeName);

    if (Bone)
    {
        Bone->Update(m_CurrentTime);
        nodeTransform = Bone->m_LocalTransform;
        OutputDebugStringA((nodeName + " local matrix (" + std::to_string(m_CurrentTime) + ") \n" + MatrixToString(nodeTransform)).c_str());
    }
    //else OutputDebugString((std::string("Did not find:") + nodeName + std::string("\n")).c_str());

    float4x4 globalTransformation = nodeTransform * parentTransform;

    auto boneInfoMap = m_CurrentAnimation->m_BoneInfoMap;
    if (boneInfoMap.find(nodeName) != boneInfoMap.end())
    {
        int index = boneInfoMap[nodeName].id;
        float4x4 offset = boneInfoMap[nodeName].transform_matrix;
        m_FinalBoneMatrices[index] = XMMatrixTranspose(offset * globalTransformation);
        //OutputDebugStringA(MatrixToString(offset).c_str());
        //OutputDebugStringA("\n");
    //    OutputDebugString((std::string("Did find*: ") + nodeName + std::string(" at index: ") + std::to_string(index) + std::string("\n")).c_str());
    }
    //else OutputDebugString((std::string("Did not find*:") + nodeName + std::string("\n")).c_str());

    for (int i = 0; i < node->childrenCount; i++)
        CalculateBoneTransform(&node->children[i], globalTransformation);
}