#include "ViewDependentStructures.fxh"

struct PixelShaderInput
{
    float4 Pos : SV_POSITION;
    float4 Color : COLOR;
};

ConstantBuffer<S_Constants> constants               : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer       : register(t0, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]        : register(t0, space1);
StructuredBuffer<S_Vertex> verticesBuffers[]        : register(t0, space2);
StructuredBuffer<uint> vertexIndicesBuffers[]       : register(t0, space3);
StructuredBuffer<uint> primitiveIndicesBuffers[]    : register(t0, space4);


uint SampleTriangleBufferAsCharArray(uint i, uint mesh_id)
{
    uint array_pos = i / 4;
    uint inner_pos = i % 4;
    uint ret = primitiveIndicesBuffers[mesh_id][array_pos];
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
          in payload S_Payload gs_Payload,
          out indices uint3 tris[MAX_MESHLET_PRIMITIVE_COUNT],
          out vertices PixelShaderInput verts[MAX_MESHLET_VERTEX_COUNT])
{
   
    
    S_PayloadEntry payload_task = gs_Payload.tasks[gid];
    S_SceneObject scene_object = objectsBuffer[payload_task.object_id];
    S_Meshlet meshlet = meshletBuffers[scene_object.mesh_id][payload_task.meshlet_id];
    
    
    SetMeshOutputCounts(meshlet.vertex_count, meshlet.triangle_count);
    
    
    // process meshlet vertices
    const uint vertexLoops = (MAX_MESHLET_VERTEX_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    for (uint v_loop = 0; v_loop < vertexLoops; v_loop++)
    {
        uint v = I + v_loop * GROUP_SIZE;
        v = min(v, meshlet.vertex_count - 1);
        
        int vertexIndex = vertexIndicesBuffers[scene_object.mesh_id][meshlet.vertex_offset + v];
        S_Vertex vertex = verticesBuffers[scene_object.mesh_id][vertexIndex];
        
        verts[v].Pos = mul(mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix), constants.ViewProjMat);
        
        float4 normal = mul(float4(vertex.normal.xyz, 0), scene_object.object_matrix);
        float brightness = clamp(clamp(dot(normalize(float3(1, 1, 1)), normal.xyz), 0, 1) + clamp(dot(normalize(float3(-2, 1, -1)), normal.xyz), 0, 0.6), 0.05, 1);
        if (constants.shadingSelection == DEFAULT_SHADING)
        {
            verts[v].Color = brightness;
        }
        else if (constants.shadingSelection == DEBUG_MESHLET_SHADING)
        {
            verts[v].Color = Rainbow(Random(gid)) * brightness;
        }
        else if (constants.shadingSelection == DEBUG_LOD_SHADING)
        {
            verts[v].Color = Rainbow(payload_task.lod_tree_depth / 6.0) * brightness;
        }
        else if (constants.shadingSelection == DEBUG_WORLD_POS)
        {
            verts[v].Color = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
        }
        else if (constants.shadingSelection == DEBUG_MESHLET_GROUP)
        {
            verts[v].Color = Rainbow(Random(payload_task.lod_morphing)) * brightness;
        }
        else if (constants.shadingSelection == DEBUG_VERTICES)
        {
            verts[v].Color = Rainbow(Random(v));
        }
        else
        {
            // should not be possible to reach??
            verts[v].Color = float4(1, 1, 0, 1);
        }
    }
    
    // process triangles
    const uint primLoops = (MAX_MESHLET_PRIMITIVE_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    for (uint p_loop = 0; p_loop < primLoops; p_loop++)
    {
        uint p = I + p_loop * GROUP_SIZE;
        p = min(p, meshlet.triangle_count - 1);

        tris[p] = uint3(SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 0, scene_object.mesh_id),
                        SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 1, scene_object.mesh_id),
                        SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 2, scene_object.mesh_id));
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    SetMeshOutputCounts(3, 1);
    verts[0].Pos = float4(-1, -1, 0, 0);
    verts[0].Color = float4(1, 0, 0, 1);
    verts[1].Pos = float4(1, -1, 0, 0);
    verts[1].Color = float4(1, 0, 0, 1);
    verts[2].Pos = float4(-1, 1, 0, 0);
    verts[2].Color = float4(1, 0, 0, 1);
    tris[0] = uint3(0, 2, 1);
    */
    
    /*
    SetMeshOutputCounts(3, 1);
    
    verts[0].Pos = mul(float4(-1, -1, 0, 1), constants.ViewProjMat);
    verts[1].Pos = mul(float4(1, -1, 0, 1), constants.ViewProjMat);
    verts[2].Pos = mul(float4(-1, 1, 0, 1), constants.ViewProjMat);
    
    verts[0].Color = float4(1, 0, 0, 1);
    verts[1].Color = float4(1, 0, 0, 1);
    verts[2].Color = float4(1, 0, 0, 1);
    
    tris[0] = uint3(0, 2, 1);
    */
    
   
    
}