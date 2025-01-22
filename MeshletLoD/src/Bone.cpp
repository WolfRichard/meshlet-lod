#include "Bone.h"
#include "Windows.h"

using namespace DirectX;

Bone::Bone(const std::string& name, int ID, const aiNodeAnim* channel)
    :
    m_Name(name),
    m_ID(ID)
{
    m_LocalTransform = DirectX::XMMatrixIdentity();

    m_NumPositions = channel->mNumPositionKeys;
    for (int positionIndex = 0; positionIndex < m_NumPositions; ++positionIndex)
    {
        aiVector3D aiPosition = channel->mPositionKeys[positionIndex].mValue;
        float timeStamp = (float)channel->mPositionKeys[positionIndex].mTime;
        KeyPosition data;
        data.position = float3(aiPosition.x, aiPosition.y, aiPosition.z);
        data.timeStamp = timeStamp;
        m_Positions.push_back(data);
    }

    m_NumRotations = channel->mNumRotationKeys;
    for (int rotationIndex = 0; rotationIndex < m_NumRotations; ++rotationIndex)
    {
        aiQuaternion aiOrientation = channel->mRotationKeys[rotationIndex].mValue;
        float timeStamp = (float)channel->mRotationKeys[rotationIndex].mTime;
        KeyRotation data;
        data.orientation = DirectX::XMVectorSet(aiOrientation.x, aiOrientation.y, aiOrientation.z, aiOrientation.w);
        data.timeStamp = timeStamp;
        m_Rotations.push_back(data);
    }

    m_NumScalings = channel->mNumScalingKeys;
    for (int keyIndex = 0; keyIndex < m_NumScalings; ++keyIndex)
    {
        aiVector3D scale = channel->mScalingKeys[keyIndex].mValue;
        float timeStamp = (float)channel->mScalingKeys[keyIndex].mTime;
        KeyScale data;
        data.scale = float3(scale.x, scale.y, scale.z);
        data.timeStamp = timeStamp;
        m_Scales.push_back(data);
    }
}

std::string MatrixToString(const DirectX::XMMATRIX& matrix) {
    std::ostringstream oss;
    const float* m = reinterpret_cast<const float*>(&matrix);  // Access matrix elements
    for (int row = 0; row < 4; ++row) {
        oss << "[ ";
        for (int col = 0; col < 4; ++col) {
            oss << m[row * 4 + col] << " ";  // Row-major order
        }
        oss << "]\n";
    }
    return oss.str();
}


void Bone::Update(float animationTime)
{
    float4x4 translation = InterpolatePosition(animationTime);
    float4x4 rotation = InterpolateRotation(animationTime);
    float4x4 scale = InterpolateScaling(animationTime);
    m_LocalTransform = scale * rotation * translation;//translation * rotation * scale;
    //OutputDebugString((MatrixToString(DirectX::XMMatrixTranslation(1, 2, 3) * DirectX::XMMatrixTranslation(1, 2, 3)) + std::string("\n\n")).c_str());
}

int Bone::GetPositionIndex(float animationTime)
{
    for (int index = 0; index < m_NumPositions - 1; ++index)
    {
        if (animationTime < m_Positions[index + 1].timeStamp)
            return index;
    }
    assert(0);
    return -1;
}

int Bone::GetRotationIndex(float animationTime)
{
    for (int index = 0; index < m_NumRotations - 1; ++index)
    {
        if (animationTime < m_Rotations[index + 1].timeStamp)
            return index;
    }
    assert(0);
    return -1;
}

int Bone::GetScaleIndex(float animationTime)
{
    for (int index = 0; index < m_NumScalings - 1; ++index)
    {
        if (animationTime < m_Scales[index + 1].timeStamp)
            return index;
    }
    assert(0);
    return -1;
}


float Bone::GetScaleFactor(float lastTimeStamp, float nextTimeStamp, float animationTime)
{
    float scaleFactor = 0.0f;
    float midWayLength = animationTime - lastTimeStamp;
    float framesDiff = nextTimeStamp - lastTimeStamp;
    scaleFactor = midWayLength / framesDiff;
    return scaleFactor;
}


float4x4 Bone::InterpolatePosition(float animationTime)
{
    if (1 == m_NumPositions)
        return DirectX::XMMatrixTranslationFromVector(DirectX::XMLoadFloat3(&(m_Positions[0].position)));
                                             

    int p0Index = GetPositionIndex(animationTime);
    int p1Index = p0Index + 1;
    float scaleFactor = GetScaleFactor(m_Positions[p0Index].timeStamp, m_Positions[p1Index].timeStamp, animationTime);
    DirectX::XMVECTOR finalPosition = DirectX::XMVectorLerp(DirectX::XMLoadFloat3(&(m_Positions[p0Index].position)), 
                                                            DirectX::XMLoadFloat3(&(m_Positions[p1Index].position)), scaleFactor);
    
    return DirectX::XMMatrixTranslationFromVector(finalPosition);
}

float4x4 Bone::InterpolateRotation(float animationTime)
{
    if (1 == m_NumRotations)
    {
        quaternion rotation = DirectX::XMQuaternionNormalize(m_Rotations[0].orientation);
        return DirectX::XMMatrixRotationQuaternion(rotation);
    }

    int p0Index = GetRotationIndex(animationTime);
    int p1Index = p0Index + 1;
    float scaleFactor = GetScaleFactor(m_Rotations[p0Index].timeStamp,
                                       m_Rotations[p1Index].timeStamp, animationTime);
    quaternion finalRotation = DirectX::XMQuaternionSlerp(m_Rotations[p0Index].orientation,
                                                          m_Rotations[p1Index].orientation, scaleFactor);
    finalRotation = DirectX::XMQuaternionNormalize(finalRotation);


    return DirectX::XMMatrixRotationQuaternion(finalRotation);
}

float4x4 Bone::InterpolateScaling(float animationTime)
{
    if (1 == m_NumScalings)
        return DirectX::XMMatrixScalingFromVector(DirectX::XMLoadFloat3(&(m_Scales[0].scale)));

    int p0Index = GetScaleIndex(animationTime);
    int p1Index = p0Index + 1;
    float scaleFactor = GetScaleFactor(m_Scales[p0Index].timeStamp,
        m_Scales[p1Index].timeStamp, animationTime);

    float3 finalScale;
    DirectX::XMStoreFloat3(&finalScale, DirectX::XMVectorLerp(DirectX::XMLoadFloat3(&(m_Scales[p0Index].scale)),
                                                              DirectX::XMLoadFloat3(&(m_Scales[p1Index].scale)), scaleFactor));

    return DirectX::XMMatrixScalingFromVector(DirectX::XMLoadFloat3(&finalScale));
}