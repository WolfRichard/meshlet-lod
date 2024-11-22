#include "Structures.fxh"

struct PixelShaderInput
{
    float4 Pos   : SV_POSITION;
    float4 Color : COLOR;
};

cbuffer InstanceIDBuffer                         : register(b0, space0)
{
    Constants constantsBuffer;
};

cbuffer InstanceIDBuffer                         : register(b1, space0)
{
    uint object_id;
};

// Objects Buffer
StructuredBuffer<SceneObject> objectsBuffer      : register(t0, space0);

// bindless buffers
StructuredBuffer<CustomVertex> verticesBuffers[] : register(t0, space1);
StructuredBuffer<uint> indicesBuffers[]          : register(t0, space2);
StructuredBuffer<uint> trianglesBuffers[]        : register(t0, space3);


uint SampleTriangleBufferAsCharArray(uint i, uint mesh_id)
{
    uint array_pos = i / 4;
    uint inner_pos = i % 4;
    uint ret = trianglesBuffers[mesh_id][array_pos];
    ret = ret >> (inner_pos * 8);
    return ret & 0x000000ff;
}


// Hash function from H. Schechter & R. Bridson, goo.gl/RXiKaH
uint Hash(uint s)
{
    s ^= 2747636419u;
    s *= 2654435769u;
    s ^= s >> 16;
    s *= 2654435769u;
    s ^= s >> 16;
    s *= 2654435769u;
    return s;
}

float Random(uint seed)
{
    return float(Hash(seed)) / 4294967295.0; // 2^32-1
}


// generate color
float4 Rainbow(float factor)
{
    float h = factor / 1.35;
    float3 col = float3(abs(h * 6.0 - 3.0) - 1.0, 2.0 - abs(h * 6.0 - 2.0), 2.0 - abs(h * 6.0 - 4.0));
    return float4(clamp(col, float3(0.0, 0.0, 0.0), float3(1.0, 1.0, 1.0)), 1.0);
}





[shader("mesh")]
[numthreads(GROUP_SIZE, 1, 1)]
[outputtopology("triangle")]
void main(in uint I : SV_GroupIndex,
          in uint gid : SV_GroupID,
          in payload Payload meshPayload,
          out indices uint3 tris[MAX_MESHLET_PRIMITIVE_COUNT],
          out vertices PixelShaderInput verts[MAX_MESHLET_VERTEX_COUNT])
{
  
    SetMeshOutputCounts(meshPayload.vertex_count[gid], meshPayload.triangle_count[gid]);
    
    const uint vertexLoops = (MAX_MESHLET_VERTEX_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    uint mesh_id = objectsBuffer[meshPayload.object_id[gid]].mesh_id;
    
    for (uint v_loop = 0; v_loop < vertexLoops; v_loop++)
    {
        uint v = I + v_loop * GROUP_SIZE;
        v = min(v, meshPayload.vertex_count[gid] - 1);
        
        int vertexIndex = indicesBuffers[mesh_id][meshPayload.vertex_offset[gid] + v];
        CustomVertex vertex = verticesBuffers[mesh_id][vertexIndex];
        verts[v].Pos = mul(mul(float4(vertex.position.xyz, 1.0), objectsBuffer[meshPayload.object_id[gid]].object_matrix), constantsBuffer.ViewProjMat);
        //verts[v].UV = vertex.uv.xy;
        
        
        float brightness = clamp(clamp(dot(normalize(float3(2, 1, -4)), vertex.normal.xyz), 0, 1) + clamp(dot(normalize(float3(-1, -2, 4)), vertex.normal.xyz), 0, 1), 0.05, 1);
        if (constantsBuffer.BoolConstants & DEBUG_VISUALS_BIT_POS)
            verts[v].Color = Rainbow(Random(gid)) * brightness;
        else 
            verts[v].Color = brightness;

    }
    
    const uint primLoops = (MAX_MESHLET_PRIMITIVE_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    
    for (uint p_loop = 0; p_loop < primLoops; p_loop++)
    {
        uint p = I + p_loop * GROUP_SIZE;
        p = min(p, meshPayload.triangle_count[gid] - 1);

        tris[p] = uint3(SampleTriangleBufferAsCharArray(meshPayload.triangle_offset[gid] + p * 3 + 0, mesh_id),
                        SampleTriangleBufferAsCharArray(meshPayload.triangle_offset[gid] + p * 3 + 1, mesh_id),
                        SampleTriangleBufferAsCharArray(meshPayload.triangle_offset[gid] + p * 3 + 2, mesh_id));
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    SetMeshOutputCounts(3, 1);
    verts[0].Pos = float4(-1, -1, 0, 0);
    verts[0].Color = float4(1, 0, 0, 1);
    verts[1].Pos = float4(1, -1, 0, 0);
    verts[1].Color = float4(1, 0, 0, 1);
    verts[2].Pos = float4(-1, 1, 0, 0);
    verts[2].Color = float4(1, 0, 0, 1);
    tris[0] = uint3(0, 1, 2);
    */
    
    
}