#pragma once

#include "HLSLnames.h"

#include "Animation.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

namespace
{
#include "../shaders/Structures.fxh"
} // namespace

class Animator
{
public:
    Animator() = default;
    void init(Animation* animation);
    void UpdateAnimation(float dt);
    void PlayAnimation(Animation* pAnimation);
    void CalculateBoneTransform(const AssimpNodeData* node, float4x4 parentTransform);

    std::vector<float4x4> m_FinalBoneMatrices;
    float m_CurrentTime;
    Animation* m_CurrentAnimation;

private:
    void recalibrateOffsets(const AssimpNodeData* node, float4x4 parentTransform);
    
    float m_DeltaTime;
};