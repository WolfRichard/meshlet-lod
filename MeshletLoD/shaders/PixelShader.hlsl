
Texture2D heightMapTexture      : register(t1, space0);
SamplerState heightMapSampler   : register(s0, space0);

struct PixelShaderInput
{
    float4 Pos   : SV_POSITION;
    float4 Color : COLOR;
    float2 UV    : TEXCOORD;
};

float4 main( PixelShaderInput IN ) : SV_Target
{
    //return float4(heightMapTexture.SampleLevel(heightMapSampler, float2(0, 0), 0).xyz, 1);
    return IN.Color;
}
