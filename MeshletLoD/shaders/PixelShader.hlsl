#include "ViewDependentStructures.fxh"

ConstantBuffer<S_Constants> constants   : register(b0, space0);

Texture2D heightMapTexture              : register(t1, space0);
SamplerState heightMapSampler           : register(s0, space0);

struct PixelShaderInput
{
    float4 Pos   : SV_POSITION;
    float4 Color : COLOR;
    float2 UV    : TEXCOORD;
};

float4 main( PixelShaderInput IN ) : SV_Target
{
    if (constants.shadingSelection == HEIGHT_MAP_SHADING)
        return float4(heightMapTexture.SampleLevel(heightMapSampler, IN.UV, 0).xyz, 1);
    else
        return IN.Color;
}
