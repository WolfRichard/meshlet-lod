#include "Structures.fxh"

struct PixelShaderInput
{
    float4 Pos : SV_POSITION;
    float4 Color : COLOR;
    float3 WPos : TEXCOORD0;
    float3 WNormal : TEXCOORD1;
    float2 UV : TEXCOORD2;
};


// read only buffers
ConstantBuffer<S_Constants> constants   : register(b0, space0);

// displacement map 
Texture2D heightMapTexture              : register(t1, space0);
SamplerState heightMapSampler           : register(s0, space0);


// Samples the displacement map according to the provided world position and normal using triplanar mapping
// Texture sampling scale and blend ratio for tilted surfaces is adjusted by the associated constants
float4 TriplanarSample(float3 worldPos, float3 worldNormal)
{
    // Normalize and abs the normal for blending weights
    float3 blend = pow(abs(worldNormal), constants.TriPlanarBlendGrade);
    blend /= dot(blend, 1.0);

    // World-space UVs
    float2 uvX = worldPos.yz * constants.TriPlanarMappingScale;
    float2 uvY = worldPos.zx * constants.TriPlanarMappingScale;
    float2 uvZ = worldPos.xy * constants.TriPlanarMappingScale;

    // Sample from each projection axis
    float4 xTex = heightMapTexture.Sample(heightMapSampler, uvX);
    float4 yTex = heightMapTexture.Sample(heightMapSampler, uvY);
    float4 zTex = heightMapTexture.Sample(heightMapSampler, uvZ);

    // Weighted blend
    return xTex * blend.x + yTex * blend.y + zTex * blend.z;
}


// pixel shader entry point
float4 main( PixelShaderInput IN ) : SV_Target
{
    // visualize displacement map if debug shading is selected, otherwise use interpolated vertex color
    if (constants.shadingSelection == HEIGHT_MAP_SHADING)
    {
        // option to use triplanar mapping for meshes that lack proper UV-cooridinates
        if (constants.BoolConstants & TRI_PLANAR_TEXTURE_MAPPING_BIT_POS)
            return TriplanarSample(IN.WPos, IN.WNormal);
        else
            return float4(heightMapTexture.SampleLevel(heightMapSampler, float2(0.001, 0.001), 0).xyz, 1);
    }
    else
        return IN.Color;
}
