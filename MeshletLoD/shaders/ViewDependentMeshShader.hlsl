#include "ViewDependentStructures.fxh"

struct PixelShaderInput
{
    float4 Pos : SV_POSITION;
    float4 Color : COLOR;
    float2 UV : TEXCOORD;
};


ConstantBuffer<S_Constants> constants               : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer       : register(t0, space0);

Texture2D heightMapTexture                          : register(t1, space0);
SamplerState heightMapSampler                       : register(s0, space0);

RWStructuredBuffer<S_PayloadEntry> payloadBuffer    : register(u2, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]        : register(t0, space1);
StructuredBuffer<S_Vertex> verticesBuffers[]        : register(t0, space2);
StructuredBuffer<uint> vertexIndicesBuffers[]       : register(t0, space3);
StructuredBuffer<uint> primitiveIndicesBuffers[]    : register(t0, space4);
StructuredBuffer<uint> morphIndicesBuffers[]        : register(t0, space5);

groupshared uint gs_vertex_count;
groupshared uint gs_triangle_count;
groupshared uint3 gs_triangles[MAX_MESHLET_PRIMITIVE_COUNT];
groupshared S_Vertex gs_vertices[MAX_MESHLET_VERTEX_COUNT];

S_Vertex barycentricVertexInterpolation(S_Vertex a, S_Vertex b, S_Vertex c, float3 w)
{
    S_Vertex result;
    result.position = a.position * w.x + b.position * w.y + c.position * w.z;
    result.color = a.color * w.x + b.color * w.y + c.color * w.z;
    result.normal = normalize(a.normal * w.x + b.normal * w.y + c.normal * w.z);
    result.uv = a.uv * w.x + b.uv * w.y + c.uv * w.z;
    return result;
}

// level 1 tesselation results in 4 triangles and 6 vertices from a single primitive
void lvl1Tessellation(S_Vertex a, S_Vertex b, S_Vertex c)
{
    uint to = 0; // triangle offset
    InterlockedAdd(gs_triangle_count, 4, to);
    uint vo = 0; // vertex offset
    InterlockedAdd(gs_vertex_count, 6, vo);
    
    gs_vertices[vo + 0] = a;
    gs_vertices[vo + 1] = b;
    gs_vertices[vo + 2] = c;
    gs_vertices[vo + 3] = barycentricVertexInterpolation(a, b, c, float3(0.5, 0.5, 0));
    gs_vertices[vo + 4] = barycentricVertexInterpolation(a, b, c, float3(0.5, 0, 0.5));
    gs_vertices[vo + 5] = barycentricVertexInterpolation(a, b, c, float3(0, 0.5, 0.5));
    
    gs_triangles[to + 0] = uint3(vo + 0, vo + 3, vo + 4);
    gs_triangles[to + 1] = uint3(vo + 3, vo + 1, vo + 5);
    gs_triangles[to + 2] = uint3(vo + 3, vo + 5, vo + 4);
    gs_triangles[to + 3] = uint3(vo + 4, vo + 5, vo + 2);
}



// level 2 tesselation results in 16 triangles and 15 vertices from a single primitive
void lvl2Tessellation(S_Vertex a, S_Vertex b, S_Vertex c)
{
    uint to = 0; // triangle offset
    InterlockedAdd(gs_triangle_count, 16, to);
    uint vo = 0; // vertex offset
    InterlockedAdd(gs_vertex_count, 15, vo);
    
    gs_vertices[vo +  0] = a;
    gs_vertices[vo +  1] = b;
    gs_vertices[vo +  2] = c;
    gs_vertices[vo +  3] = barycentricVertexInterpolation(a, b, c, float3(0.50, 0.50, 0.00));
    gs_vertices[vo +  4] = barycentricVertexInterpolation(a, b, c, float3(0.50, 0.00, 0.50));
    gs_vertices[vo +  5] = barycentricVertexInterpolation(a, b, c, float3(0.00, 0.50, 0.50));
    gs_vertices[vo +  6] = barycentricVertexInterpolation(a, b, c, float3(0.75, 0.00, 0.25));
    gs_vertices[vo +  7] = barycentricVertexInterpolation(a, b, c, float3(0.75, 0.25, 0.00));
    gs_vertices[vo +  8] = barycentricVertexInterpolation(a, b, c, float3(0.50, 0.25, 0.25));
    gs_vertices[vo +  9] = barycentricVertexInterpolation(a, b, c, float3(0.25, 0.00, 0.75));
    gs_vertices[vo + 10] = barycentricVertexInterpolation(a, b, c, float3(0.25, 0.25, 0.50));
    gs_vertices[vo + 11] = barycentricVertexInterpolation(a, b, c, float3(0.25, 0.50, 0.25));
    gs_vertices[vo + 12] = barycentricVertexInterpolation(a, b, c, float3(0.25, 0.75, 0.00));
    gs_vertices[vo + 13] = barycentricVertexInterpolation(a, b, c, float3(0.00, 0.25, 0.75));
    gs_vertices[vo + 14] = barycentricVertexInterpolation(a, b, c, float3(0.00, 0.75, 0.25));
    
    
    gs_triangles[to +  0] = uint3( 0,  7,  6);
    gs_triangles[to +  1] = uint3( 6,  8,  4);
    gs_triangles[to +  2] = uint3( 6,  7,  8);
    gs_triangles[to +  3] = uint3( 7,  3,  8);
    gs_triangles[to +  4] = uint3( 4, 10,  9);
    gs_triangles[to +  5] = uint3( 4,  8, 10);
    gs_triangles[to +  6] = uint3( 8, 11, 10);
    gs_triangles[to +  7] = uint3( 8,  3, 11);
    gs_triangles[to +  8] = uint3( 3, 12, 11);
    gs_triangles[to +  9] = uint3( 9, 13,  2);
    gs_triangles[to + 10] = uint3( 9, 10, 13);
    gs_triangles[to + 11] = uint3(10,  5, 13);
    gs_triangles[to + 12] = uint3(10, 11,  5);
    gs_triangles[to + 13] = uint3(11, 14,  5);
    gs_triangles[to + 14] = uint3(11, 12, 14);
    gs_triangles[to + 15] = uint3(12,  1, 14);
}


/* level 3 tesselation results in 64 triangles and 45 vertices from a single primitive

             44(c)
            /   \
          ... -- ...
         /   \  /   \
        9 --  10 -- ...
       /  \ /  \  /   \   
   0(a) -- 1 -- ... -- 8(b)

*/


void lvl3Tessellation(S_Vertex a, S_Vertex b, S_Vertex c)
{
    uint to = 0; // triangle offset
    InterlockedAdd(gs_triangle_count, 64, to);
    uint vo = 0; // vertex offset
    InterlockedAdd(gs_vertex_count, 45, vo);
    
    gs_vertices[vo +  0] = a;
    gs_vertices[vo +  1] = barycentricVertexInterpolation(a, b, c, float3(0.875, 0.125, 0.000));
    gs_vertices[vo +  2] = barycentricVertexInterpolation(a, b, c, float3(0.750, 0.250, 0.000));
    gs_vertices[vo +  3] = barycentricVertexInterpolation(a, b, c, float3(0.625, 0.375, 0.000));
    gs_vertices[vo +  4] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.500, 0.000));
    gs_vertices[vo +  5] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.615, 0.000));
    gs_vertices[vo +  6] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.750, 0.000));
    gs_vertices[vo +  7] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.875, 0.000));
    gs_vertices[vo +  8] = b;
    gs_vertices[vo +  9] = barycentricVertexInterpolation(a, b, c, float3(0.875, 0.000, 0.125));
    gs_vertices[vo + 10] = barycentricVertexInterpolation(a, b, c, float3(0.750, 0.125, 0.125));
    gs_vertices[vo + 11] = barycentricVertexInterpolation(a, b, c, float3(0.625, 0.250, 0.125));
    gs_vertices[vo + 12] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.375, 0.125));
    gs_vertices[vo + 13] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.500, 0.125));
    gs_vertices[vo + 14] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.625, 0.125));
    gs_vertices[vo + 15] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.750, 0.125));
    gs_vertices[vo + 16] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.875, 0.125));
    gs_vertices[vo + 17] = barycentricVertexInterpolation(a, b, c, float3(0.750, 0.000, 0.250));
    gs_vertices[vo + 18] = barycentricVertexInterpolation(a, b, c, float3(0.625, 0.125, 0.250));
    gs_vertices[vo + 19] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.250, 0.250));
    gs_vertices[vo + 20] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.375, 0.250));
    gs_vertices[vo + 21] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.500, 0.250));
    gs_vertices[vo + 22] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.625, 0.250));
    gs_vertices[vo + 23] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.750, 0.250));
    gs_vertices[vo + 24] = barycentricVertexInterpolation(a, b, c, float3(0.625, 0.000, 0.375));
    gs_vertices[vo + 25] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.125, 0.375));
    gs_vertices[vo + 26] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.250, 0.375));
    gs_vertices[vo + 27] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.375, 0.375));
    gs_vertices[vo + 28] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.500, 0.375));
    gs_vertices[vo + 29] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.625, 0.375));
    gs_vertices[vo + 30] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.000, 0.500));
    gs_vertices[vo + 31] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.125, 0.500));
    gs_vertices[vo + 32] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.250, 0.500));
    gs_vertices[vo + 33] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.375, 0.500));
    gs_vertices[vo + 34] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.500, 0.500));
    gs_vertices[vo + 35] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.000, 0.625));
    gs_vertices[vo + 36] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.125, 0.625));
    gs_vertices[vo + 37] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.250, 0.625));
    gs_vertices[vo + 38] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.375, 0.625));
    gs_vertices[vo + 39] = barycentricVertexInterpolation(a, b, c, float3(0.250, 0.000, 0.750));
    gs_vertices[vo + 40] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.125, 0.750));
    gs_vertices[vo + 41] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.250, 0.750));
    gs_vertices[vo + 42] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.000, 0.875));
    gs_vertices[vo + 43] = barycentricVertexInterpolation(a, b, c, float3(0.125, 0.000, 0.875));
    gs_vertices[vo + 44] = c;
    
    
    gs_triangles[to +  0] = uint3(vo +  0, vo +  1, vo +  9);
    gs_triangles[to +  1] = uint3(vo +  1, vo + 10, vo +  9);
    gs_triangles[to +  2] = uint3(vo +  1, vo +  2, vo + 10);
    gs_triangles[to +  3] = uint3(vo +  2, vo + 11, vo + 10);
    gs_triangles[to +  4] = uint3(vo +  2, vo +  3, vo + 11);
    gs_triangles[to +  5] = uint3(vo +  3, vo + 12, vo + 11);
    gs_triangles[to +  6] = uint3(vo +  3, vo +  4, vo + 12);
    gs_triangles[to +  7] = uint3(vo +  4, vo + 13, vo + 12);
    gs_triangles[to +  8] = uint3(vo +  4, vo +  5, vo + 13);
    gs_triangles[to +  9] = uint3(vo +  5, vo + 14, vo + 13);
    gs_triangles[to + 10] = uint3(vo +  5, vo +  6, vo + 14);
    gs_triangles[to + 11] = uint3(vo +  6, vo + 15, vo + 14);
    gs_triangles[to + 12] = uint3(vo +  6, vo +  7, vo + 15);
    gs_triangles[to + 13] = uint3(vo +  7, vo + 16, vo + 15);
    gs_triangles[to + 14] = uint3(vo +  7, vo +  8, vo + 16);
    gs_triangles[to + 15] = uint3(vo +  9, vo + 10, vo + 17);
    gs_triangles[to + 16] = uint3(vo + 10, vo + 18, vo + 17);
    gs_triangles[to + 17] = uint3(vo + 10, vo + 11, vo + 18);
    gs_triangles[to + 18] = uint3(vo + 11, vo + 19, vo + 18);
    gs_triangles[to + 19] = uint3(vo + 11, vo + 12, vo + 19);
    gs_triangles[to + 20] = uint3(vo + 12, vo + 20, vo + 19);
    gs_triangles[to + 21] = uint3(vo + 12, vo + 13, vo + 20);
    gs_triangles[to + 22] = uint3(vo + 13, vo + 21, vo + 20);
    gs_triangles[to + 23] = uint3(vo + 13, vo + 14, vo + 21);
    gs_triangles[to + 24] = uint3(vo + 14, vo + 22, vo + 21);
    gs_triangles[to + 25] = uint3(vo + 14, vo + 15, vo + 22);
    gs_triangles[to + 26] = uint3(vo + 15, vo + 23, vo + 22);
    gs_triangles[to + 27] = uint3(vo + 15, vo + 16, vo + 23);
    gs_triangles[to + 28] = uint3(vo + 17, vo + 18, vo + 24);
    gs_triangles[to + 29] = uint3(vo + 18, vo + 25, vo + 24);
    gs_triangles[to + 30] = uint3(vo + 18, vo + 19, vo + 25);
    gs_triangles[to + 31] = uint3(vo + 19, vo + 26, vo + 25);
    gs_triangles[to + 32] = uint3(vo + 19, vo + 20, vo + 26);
    gs_triangles[to + 33] = uint3(vo + 20, vo + 27, vo + 26);
    gs_triangles[to + 34] = uint3(vo + 20, vo + 21, vo + 27);
    gs_triangles[to + 35] = uint3(vo + 21, vo + 28, vo + 27);
    gs_triangles[to + 36] = uint3(vo + 21, vo + 22, vo + 28);
    gs_triangles[to + 37] = uint3(vo + 22, vo + 29, vo + 28);
    gs_triangles[to + 38] = uint3(vo + 22, vo + 23, vo + 29);
    gs_triangles[to + 39] = uint3(vo + 24, vo + 25, vo + 30);
    gs_triangles[to + 40] = uint3(vo + 25, vo + 31, vo + 30);
    gs_triangles[to + 41] = uint3(vo + 25, vo + 26, vo + 31);
    gs_triangles[to + 42] = uint3(vo + 26, vo + 32, vo + 31);
    gs_triangles[to + 43] = uint3(vo + 26, vo + 27, vo + 32);
    gs_triangles[to + 44] = uint3(vo + 27, vo + 33, vo + 32);
    gs_triangles[to + 45] = uint3(vo + 27, vo + 28, vo + 33);
    gs_triangles[to + 46] = uint3(vo + 28, vo + 34, vo + 33);
    gs_triangles[to + 47] = uint3(vo + 28, vo + 29, vo + 34);
    gs_triangles[to + 48] = uint3(vo + 30, vo + 31, vo + 35);
    gs_triangles[to + 49] = uint3(vo + 31, vo + 36, vo + 35);
    gs_triangles[to + 50] = uint3(vo + 31, vo + 32, vo + 36);
    gs_triangles[to + 51] = uint3(vo + 32, vo + 37, vo + 36);
    gs_triangles[to + 52] = uint3(vo + 32, vo + 33, vo + 37);
    gs_triangles[to + 53] = uint3(vo + 33, vo + 38, vo + 37);
    gs_triangles[to + 54] = uint3(vo + 33, vo + 34, vo + 38);
    gs_triangles[to + 55] = uint3(vo + 35, vo + 36, vo + 39);
    gs_triangles[to + 56] = uint3(vo + 36, vo + 40, vo + 39);
    gs_triangles[to + 57] = uint3(vo + 36, vo + 37, vo + 40);
    gs_triangles[to + 58] = uint3(vo + 37, vo + 41, vo + 40);
    gs_triangles[to + 59] = uint3(vo + 37, vo + 38, vo + 41);
    gs_triangles[to + 60] = uint3(vo + 39, vo + 40, vo + 42);
    gs_triangles[to + 61] = uint3(vo + 40, vo + 43, vo + 42);
    gs_triangles[to + 62] = uint3(vo + 40, vo + 41, vo + 43);
    gs_triangles[to + 63] = uint3(vo + 42, vo + 43, vo + 44);
}


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

float getExpectedLoDLevel(float4 position)// position must be in world space!
{
    float cam_dist = distance(constants.CameraWorldPos, position.xyz);
    return max(log2(cam_dist / constants.LoD_Scale), 0);
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
   
    
    S_PayloadEntry payload_task = payloadBuffer[gid + gs_Payload.global_payload_offset]; //    gs_Payload.tasks[gid];

    S_SceneObject scene_object = objectsBuffer[payload_task.object_id];
    S_Meshlet meshlet = meshletBuffers[scene_object.mesh_id][payload_task.meshlet_id];
    
    
    /*
    
    // optional tessellation
    bool do_tessellation = ((constants.BoolConstants & TRESSELLATION_BIT_POS) && payload_task.tessellation_grade);
    if (do_tessellation)
    {
        // Reset the group shared counter variables from the first thread in the group
        if (I == 0)
        {
            gs_MeshletCount = 0;
        }
        GroupMemoryBarrierWithGroupSync();
        
        // tessellate each individual triangle
        for (uint p = payload_task.tessellation_triangle_offset + I; p < payload_task.tessellation_triangle_count + payload_task.tessellation_triangle_offset; p += GROUP_SIZE)
        {
            S_Vertex a = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 0, scene_object.mesh_id)]];
            S_Vertex b = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 1, scene_object.mesh_id)]];
            S_Vertex c = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 2, scene_object.mesh_id)]];
        
            if (payload_task.tessellation_grade == 3)
                lvl3Tessellation(a, b, c);
            else if ( payload_task.tessellation_grade == 2)
                lvl2Tessellation(a, b, c);
            else
                lvl1Tessellation(a, b, c);
        }
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    */
    
    
    SetMeshOutputCounts(meshlet.vertex_count, meshlet.triangle_count);
    
    // process meshlet vertices
    const uint vertexLoops = (MAX_MESHLET_VERTEX_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    for (uint v_loop = 0; v_loop < vertexLoops; v_loop++)
    {
        uint v = I + v_loop * GROUP_SIZE;
        v = min(v, meshlet.vertex_count - 1);
        
        int vertexIndex = vertexIndicesBuffers[scene_object.mesh_id][meshlet.vertex_offset + v];
        S_Vertex vertex = verticesBuffers[scene_object.mesh_id][vertexIndex];
        
        float4 original_world_pos = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
        
        float expected_LoD = getExpectedLoDLevel(original_world_pos);
        float lerp_value = clamp(expected_LoD - meshlet.discrete_level_of_detail, -100, 100);
        int morphTargetIndex = morphIndicesBuffers[scene_object.mesh_id][meshlet.vertex_offset + v];
        
        if ((constants.BoolConstants & GEO_MORPHING_BIT_POS) && !(constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS))
        {
            while (lerp_value > 1)
            {
                vertexIndex = vertexIndicesBuffers[scene_object.mesh_id][morphTargetIndex];
                morphTargetIndex = morphIndicesBuffers[scene_object.mesh_id][morphTargetIndex];
                lerp_value -= 1;
                
            }
            vertex = verticesBuffers[scene_object.mesh_id][vertexIndex];
            //int morphTargetIndex = morphIndicesBuffers[scene_object.mesh_id][meshlet.vertex_offset + v];
            int morphTargetVertexIndex = vertexIndicesBuffers[scene_object.mesh_id][morphTargetIndex];
            S_Vertex morphTargetVertex = verticesBuffers[scene_object.mesh_id][morphTargetVertexIndex];
              
            vertex.position = lerp(vertex.position, morphTargetVertex.position, lerp_value);
            vertex.normal   = lerp(vertex.normal,   morphTargetVertex.normal,   lerp_value);
            vertex.color    = lerp(vertex.color,    morphTargetVertex.color,    lerp_value);
            vertex.uv       = lerp(vertex.uv,       morphTargetVertex.uv,       lerp_value);
            
            vertex.normal   = normalize(vertex.normal);
        }
        
        
        verts[v].Pos = mul(mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix), constants.ViewProjMat);
        verts[v].UV = vertex.uv.xy;
        
        float4 normal = mul(float4(vertex.normal.xyz, 0), scene_object.object_matrix);
        float brightness = clamp(clamp(dot(normalize(float3(1, 1, 1)), normal.xyz), 0, 1) + clamp(dot(normalize(float3(-2, 1, -1)), normal.xyz), 0, 0.6), 0.05, 1);
        if (constants.shadingSelection == DEFAULT_SHADING)
        {
            verts[v].Color = brightness;
        }
        else if (constants.shadingSelection == MESHLETS_SHADING)
        {
            verts[v].Color = Rainbow(Random(payload_task.meshlet_id)) * brightness;
        }
        else if (constants.shadingSelection == LOD_SHADING)
        {
            if ((constants.BoolConstants & GEO_MORPHING_BIT_POS) && !(constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS))
                verts[v].Color = Rainbow(expected_LoD / 5.0);
            else
                verts[v].Color = Rainbow(meshlet.discrete_level_of_detail / 5.0) * brightness;
            
        }
        else if (constants.shadingSelection == MESHLET_GROUP_SHADING)
        {
            verts[v].Color = Rainbow(Random(meshlet.group_id)) * brightness;
        }
        else if (constants.shadingSelection == TESSELLATION_LEVEL_SHADING)
        {
            verts[v].Color = Rainbow(payload_task.tessellation_grade / 2.75) * brightness;
        }
        else
        {
            // should not be shown. resulting color is probably going to be overwritten by pixel shader
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