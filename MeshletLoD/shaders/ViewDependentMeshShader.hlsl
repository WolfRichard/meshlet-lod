#include "ViewDependentStructures.fxh"

struct PixelShaderInput
{
    float4 Pos      : SV_POSITION;
    float4 Color    : COLOR;
    float3 WPos     : TEXCOORD0;
    float3 WNormal  : TEXCOORD1;
    float2 UV       : TEXCOORD2;
};

// read only buffers
ConstantBuffer<S_Constants> constants               : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer       : register(t0, space0);

// displacement map 
Texture2D heightMapTexture                          : register(t1, space0);
SamplerState heightMapSampler                       : register(s0, space0);

// read-write buffer
RWStructuredBuffer<S_PayloadEntry> payloadBuffer    : register(u2, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]        : register(t0, space1);
StructuredBuffer<S_Vertex> verticesBuffers[]        : register(t0, space2);
StructuredBuffer<uint> vertexIndicesBuffers[]       : register(t0, space3);
StructuredBuffer<uint> primitiveIndicesBuffers[]    : register(t0, space4);
StructuredBuffer<uint> morphIndicesBuffers[]        : register(t0, space5);

// group shared variables
groupshared uint gs_vertex_count;
groupshared uint gs_triangle_count;
groupshared uint3 gs_triangles[MAX_MESHLET_PRIMITIVE_COUNT];
groupshared S_Vertex gs_vertices[MAX_MESHLET_VERTEX_COUNT];

groupshared uint gs_tessellation_levels[MAX_MESHLET_VERTEX_COUNT];          // per vertex tessellation level (even lvl3 tessellation produces lvl1 and lvl2 tessellation vertices) 
groupshared uint gs_tessellation_morph_targets[MAX_MESHLET_VERTEX_COUNT];   // indices into what tessellation vertex the current vertex should be interpolated towards when applying geo-morphing.  for tessellation lvl0 vertices (base vertices of the triangle, the morph target index is the original morph target index of the morph-target global buffer)


// helper function that linearly interpolates all components of a S_Vertex
S_Vertex linearVertexInterpolation(S_Vertex a, S_Vertex b, float w)
{
    S_Vertex result;
    result.position = lerp(a.position, b.position, w);
    result.color    = lerp(a.color   , b.color   , w);
    result.normal   = lerp(a.normal  , b.normal  , w);
    result.uv       = lerp(a.uv      , b.uv      , w);
    
    result.normal   = normalize(result.normal);
    return result;
}


// helper function that uses barycentric coordinates to interpolate all components of a S_Vertex
S_Vertex barycentricVertexInterpolation(S_Vertex a, S_Vertex b, S_Vertex c, float3 w)
{
    S_Vertex result;
    
    result.position = a.position * w.x + b.position * w.y + c.position * w.z;
    result.color    = a.color    * w.x + b.color    * w.y + c.color    * w.z;
    result.normal   = a.normal   * w.x + b.normal   * w.y + c.normal   * w.z;
    result.uv       = a.uv       * w.x + b.uv       * w.y + c.uv       * w.z;
    
    result.normal   = normalize(result.normal);
    return result;
}



// level 1 tesselation results in 4 triangles and 6 vertices from a single primitive
void lvl1Tessellation(uint3 vertex_indices, uint mesh_id)
{
    uint index_a = vertexIndicesBuffers[mesh_id][vertex_indices.x];
    uint index_b = vertexIndicesBuffers[mesh_id][vertex_indices.y];
    uint index_c = vertexIndicesBuffers[mesh_id][vertex_indices.z];
    
    S_Vertex a = verticesBuffers[mesh_id][index_a];
    S_Vertex b = verticesBuffers[mesh_id][index_b];
    S_Vertex c = verticesBuffers[mesh_id][index_c];
    
    uint morph_target_a = morphIndicesBuffers[mesh_id][vertex_indices.x];
    uint morph_target_b = morphIndicesBuffers[mesh_id][vertex_indices.y];
    uint morph_target_c = morphIndicesBuffers[mesh_id][vertex_indices.z];
    
    
    uint to = 0; // triangle offset
    InterlockedAdd(gs_triangle_count, 4, to);
    uint vo = 0; // vertex offset
    InterlockedAdd(gs_vertex_count, 6, vo);
    
    // generate tesselated vertices
    gs_vertices[vo + 0] = a;
    gs_vertices[vo + 1] = b;
    gs_vertices[vo + 2] = c;
    gs_vertices[vo + 3] = barycentricVertexInterpolation(a, b, c, float3(0.5, 0.5, 0));
    gs_vertices[vo + 4] = barycentricVertexInterpolation(a, b, c, float3(0.5, 0, 0.5));
    gs_vertices[vo + 5] = barycentricVertexInterpolation(a, b, c, float3(0, 0.5, 0.5));
    
    //generate connectivity between newly produced vertices
    gs_triangles[to + 0] = uint3(vo + 0, vo + 3, vo + 4);
    gs_triangles[to + 1] = uint3(vo + 3, vo + 1, vo + 5);
    gs_triangles[to + 2] = uint3(vo + 3, vo + 5, vo + 4);
    gs_triangles[to + 3] = uint3(vo + 4, vo + 5, vo + 2);
    
    
    // generation of Geo-Morph information
    
    // check what winding order edge tessellation should be morphed to, to avoid that opposing T-junktions drift appart they must decide to morph into the same vertices
    // always morph to the higher index, default assumption is that ordered by hash A > B > C > A. (this can/should never be true, and just represents their original winding order)
    // inverse edge morphing for specific edges if indices order deveates from default assumption
    bool i_ab = index_a < index_b;
    bool i_bc = index_b < index_c;
    bool i_ca = index_c < index_a;
    
    gs_tessellation_levels[vo + 0] = 0;
    gs_tessellation_levels[vo + 1] = 0;
    gs_tessellation_levels[vo + 2] = 0;
    gs_tessellation_levels[vo + 3] = 1;
    gs_tessellation_levels[vo + 4] = 1;
    gs_tessellation_levels[vo + 5] = 1;
    
    gs_tessellation_morph_targets[vo + 0] = morph_target_a;
    gs_tessellation_morph_targets[vo + 1] = morph_target_b;
    gs_tessellation_morph_targets[vo + 2] = morph_target_c;
    gs_tessellation_morph_targets[vo + 3] = vo + (i_ab ? 1 : 0);
    gs_tessellation_morph_targets[vo + 4] = vo + (i_ca ? 0 : 2);
    gs_tessellation_morph_targets[vo + 5] = vo + (i_bc ? 2 : 1);
}



// level 2 tesselation results in 16 triangles and 15 vertices from a single primitive
void lvl2Tessellation(uint3 vertex_indices, uint mesh_id)
{
    uint index_a = vertexIndicesBuffers[mesh_id][vertex_indices.x];
    uint index_b = vertexIndicesBuffers[mesh_id][vertex_indices.y];
    uint index_c = vertexIndicesBuffers[mesh_id][vertex_indices.z];
    
    S_Vertex a = verticesBuffers[mesh_id][index_a];
    S_Vertex b = verticesBuffers[mesh_id][index_b];
    S_Vertex c = verticesBuffers[mesh_id][index_c];
    
    uint morph_target_a = morphIndicesBuffers[mesh_id][vertex_indices.x];
    uint morph_target_b = morphIndicesBuffers[mesh_id][vertex_indices.y];
    uint morph_target_c = morphIndicesBuffers[mesh_id][vertex_indices.z];
    
    
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
    
    
    gs_triangles[to +  0] = uint3(vo +  0, vo +  7, vo +  6);
    gs_triangles[to +  1] = uint3(vo +  6, vo +  8, vo +  4);
    gs_triangles[to +  2] = uint3(vo +  6, vo +  7, vo +  8);
    gs_triangles[to +  3] = uint3(vo +  7, vo +  3, vo +  8);
    gs_triangles[to +  4] = uint3(vo +  4, vo + 10, vo +  9);
    gs_triangles[to +  5] = uint3(vo +  4, vo +  8, vo + 10);
    gs_triangles[to +  6] = uint3(vo +  8, vo + 11, vo + 10);
    gs_triangles[to +  7] = uint3(vo +  8, vo +  3, vo + 11);
    gs_triangles[to +  8] = uint3(vo +  3, vo + 12, vo + 11);
    gs_triangles[to +  9] = uint3(vo +  9, vo + 13, vo +  2);
    gs_triangles[to + 10] = uint3(vo +  9, vo + 10, vo + 13);
    gs_triangles[to + 11] = uint3(vo + 10, vo +  5, vo + 13);
    gs_triangles[to + 12] = uint3(vo + 10, vo + 11, vo +  5);
    gs_triangles[to + 13] = uint3(vo + 11, vo + 14, vo +  5);
    gs_triangles[to + 14] = uint3(vo + 11, vo + 12, vo + 14);
    gs_triangles[to + 15] = uint3(vo + 12, vo +  1, vo + 14);
    
    
    // generation of Geo-Morph informatin
    bool i_ab = index_a < index_b;
    bool i_bc = index_b < index_c;
    bool i_ca = index_c < index_a;
    
    gs_tessellation_levels[vo +  0] = 0;
    gs_tessellation_levels[vo +  1] = 0;
    gs_tessellation_levels[vo +  2] = 0;
    gs_tessellation_levels[vo +  3] = 1;
    gs_tessellation_levels[vo +  4] = 1;
    gs_tessellation_levels[vo +  5] = 1;
    gs_tessellation_levels[vo +  6] = 2;
    gs_tessellation_levels[vo +  7] = 2;
    gs_tessellation_levels[vo +  8] = 2;
    gs_tessellation_levels[vo +  9] = 2;
    gs_tessellation_levels[vo + 10] = 2;
    gs_tessellation_levels[vo + 11] = 2;
    gs_tessellation_levels[vo + 12] = 2;
    gs_tessellation_levels[vo + 13] = 2;
    gs_tessellation_levels[vo + 14] = 2;
    
    
    gs_tessellation_morph_targets[vo +  0] = morph_target_a;
    gs_tessellation_morph_targets[vo +  1] = morph_target_b;
    gs_tessellation_morph_targets[vo +  2] = morph_target_c;
    gs_tessellation_morph_targets[vo +  3] = vo + (i_ab ? 1 : 0);
    gs_tessellation_morph_targets[vo +  4] = vo + (i_ca ? 0 : 2);
    gs_tessellation_morph_targets[vo +  5] = vo + (i_bc ? 2 : 1);
    gs_tessellation_morph_targets[vo +  6] = vo + 4;
    gs_tessellation_morph_targets[vo +  7] = vo + 3;
    gs_tessellation_morph_targets[vo +  8] = vo + 3;
    gs_tessellation_morph_targets[vo +  9] = vo + 4;
    gs_tessellation_morph_targets[vo + 10] = vo + 4;
    gs_tessellation_morph_targets[vo + 11] = vo + 5;
    gs_tessellation_morph_targets[vo + 12] = vo + 3;
    gs_tessellation_morph_targets[vo + 13] = vo + 5;
    gs_tessellation_morph_targets[vo + 14] = vo + 5; 
}


/* level 3 tesselation results in 64 triangles and 45 vertices from a single primitive
   TODO: this is quite a heavy load for a single thread, this should be distributed among all threads in the future

             44(c)
            /   \
          ... -- ...
         /   \  /   \
        9 --  10 -- ...
       /  \ /  \  /   \   
   0(a) -- 1 -- ... -- 8(b)

*/


void lvl3Tessellation(uint3 vertex_indices, uint mesh_id)
{
    uint index_a = vertexIndicesBuffers[mesh_id][vertex_indices.x];
    uint index_b = vertexIndicesBuffers[mesh_id][vertex_indices.y];
    uint index_c = vertexIndicesBuffers[mesh_id][vertex_indices.z];
    
    S_Vertex a = verticesBuffers[mesh_id][index_a];
    S_Vertex b = verticesBuffers[mesh_id][index_b];
    S_Vertex c = verticesBuffers[mesh_id][index_c];
    
    uint morph_target_a = morphIndicesBuffers[mesh_id][vertex_indices.x];
    uint morph_target_b = morphIndicesBuffers[mesh_id][vertex_indices.y];
    uint morph_target_c = morphIndicesBuffers[mesh_id][vertex_indices.z];
    
    
    uint to = 0; // triangle offset
    InterlockedAdd(gs_triangle_count, 64, to);
    uint vo = 0; // vertex offset
    InterlockedAdd(gs_vertex_count, 45, vo);
    
    gs_vertices[vo +  0] = a;
    gs_vertices[vo +  1] = barycentricVertexInterpolation(a, b, c, float3(0.875, 0.125, 0.000));
    gs_vertices[vo +  2] = barycentricVertexInterpolation(a, b, c, float3(0.750, 0.250, 0.000));
    gs_vertices[vo +  3] = barycentricVertexInterpolation(a, b, c, float3(0.625, 0.375, 0.000));
    gs_vertices[vo +  4] = barycentricVertexInterpolation(a, b, c, float3(0.500, 0.500, 0.000));
    gs_vertices[vo +  5] = barycentricVertexInterpolation(a, b, c, float3(0.375, 0.625, 0.000));
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
    gs_vertices[vo + 43] = barycentricVertexInterpolation(a, b, c, float3(0.000, 0.125, 0.875));
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
    
    // generation of Geo-Morph informatin
    bool i_ab = index_a < index_b;
    bool i_bc = index_b < index_c;
    bool i_ca = index_c < index_a;
    
    gs_tessellation_levels[vo +  0] = 0;
    gs_tessellation_levels[vo +  1] = 3;
    gs_tessellation_levels[vo +  2] = 2;
    gs_tessellation_levels[vo +  3] = 3;
    gs_tessellation_levels[vo +  4] = 1;
    gs_tessellation_levels[vo +  5] = 3;
    gs_tessellation_levels[vo +  6] = 2;
    gs_tessellation_levels[vo +  7] = 3;
    gs_tessellation_levels[vo +  8] = 0;
    gs_tessellation_levels[vo +  9] = 3;
    gs_tessellation_levels[vo + 10] = 3;
    gs_tessellation_levels[vo + 11] = 3;
    gs_tessellation_levels[vo + 12] = 3;
    gs_tessellation_levels[vo + 13] = 3;
    gs_tessellation_levels[vo + 14] = 3;
    gs_tessellation_levels[vo + 15] = 3;
    gs_tessellation_levels[vo + 16] = 3;
    gs_tessellation_levels[vo + 17] = 2;
    gs_tessellation_levels[vo + 18] = 3;
    gs_tessellation_levels[vo + 19] = 2;
    gs_tessellation_levels[vo + 20] = 3;
    gs_tessellation_levels[vo + 21] = 2;
    gs_tessellation_levels[vo + 22] = 3;
    gs_tessellation_levels[vo + 23] = 2;
    gs_tessellation_levels[vo + 24] = 3;
    gs_tessellation_levels[vo + 25] = 3;
    gs_tessellation_levels[vo + 26] = 3;
    gs_tessellation_levels[vo + 27] = 3;
    gs_tessellation_levels[vo + 28] = 3;
    gs_tessellation_levels[vo + 29] = 3;
    gs_tessellation_levels[vo + 30] = 1;
    gs_tessellation_levels[vo + 31] = 3;
    gs_tessellation_levels[vo + 32] = 2;
    gs_tessellation_levels[vo + 33] = 3;
    gs_tessellation_levels[vo + 34] = 1;
    gs_tessellation_levels[vo + 35] = 3;
    gs_tessellation_levels[vo + 36] = 3;
    gs_tessellation_levels[vo + 37] = 3;
    gs_tessellation_levels[vo + 38] = 3;
    gs_tessellation_levels[vo + 39] = 2;
    gs_tessellation_levels[vo + 40] = 3;
    gs_tessellation_levels[vo + 41] = 2;
    gs_tessellation_levels[vo + 42] = 3;
    gs_tessellation_levels[vo + 43] = 3;
    gs_tessellation_levels[vo + 44] = 0;
    
    gs_tessellation_morph_targets[vo +  0] = morph_target_a;
    gs_tessellation_morph_targets[vo +  1] = vo + 2;
    gs_tessellation_morph_targets[vo +  2] = vo + 4;
    gs_tessellation_morph_targets[vo +  3] = vo + 2;
    gs_tessellation_morph_targets[vo +  4] = vo + (i_ab ? 8 : 0);
    gs_tessellation_morph_targets[vo +  5] = vo + 6;
    gs_tessellation_morph_targets[vo +  6] = vo + 4;
    gs_tessellation_morph_targets[vo +  7] = vo + 6;
    gs_tessellation_morph_targets[vo +  8] = morph_target_b;
    gs_tessellation_morph_targets[vo +  9] = vo + 17;
    gs_tessellation_morph_targets[vo + 10] = vo + 17;
    gs_tessellation_morph_targets[vo + 11] = vo + 2;
    gs_tessellation_morph_targets[vo + 12] = vo + 19;
    gs_tessellation_morph_targets[vo + 13] = vo + 21;
    gs_tessellation_morph_targets[vo + 14] = vo + 21;
    gs_tessellation_morph_targets[vo + 15] = vo + 6;
    gs_tessellation_morph_targets[vo + 16] = vo + 23;
    gs_tessellation_morph_targets[vo + 17] = vo + 30;
    gs_tessellation_morph_targets[vo + 18] = vo + 19;
    gs_tessellation_morph_targets[vo + 19] = vo + 4;
    gs_tessellation_morph_targets[vo + 20] = vo + 21;
    gs_tessellation_morph_targets[vo + 21] = vo + 34;
    gs_tessellation_morph_targets[vo + 22] = vo + 23;
    gs_tessellation_morph_targets[vo + 23] = vo + 34;
    gs_tessellation_morph_targets[vo + 24] = vo + 17;
    gs_tessellation_morph_targets[vo + 25] = vo + 19;
    gs_tessellation_morph_targets[vo + 26] = vo + 19;
    gs_tessellation_morph_targets[vo + 27] = vo + 32; 
    gs_tessellation_morph_targets[vo + 28] = vo + 21;
    gs_tessellation_morph_targets[vo + 29] = vo + 23;
    gs_tessellation_morph_targets[vo + 30] = vo + (i_ca ? 0 : 44);
    gs_tessellation_morph_targets[vo + 31] = vo + 32;
    gs_tessellation_morph_targets[vo + 32] = vo + 30;
    gs_tessellation_morph_targets[vo + 33] = vo + 32;
    gs_tessellation_morph_targets[vo + 34] = vo + (i_bc ? 44 : 8);
    gs_tessellation_morph_targets[vo + 35] = vo + 39;
    gs_tessellation_morph_targets[vo + 36] = vo + 39;
    gs_tessellation_morph_targets[vo + 37] = vo + 32;
    gs_tessellation_morph_targets[vo + 38] = vo + 41;
    gs_tessellation_morph_targets[vo + 39] = vo + 30;
    gs_tessellation_morph_targets[vo + 40] = vo + 41;
    gs_tessellation_morph_targets[vo + 41] = vo + 34;
    gs_tessellation_morph_targets[vo + 42] = vo + 39;
    gs_tessellation_morph_targets[vo + 43] = vo + 41;
    gs_tessellation_morph_targets[vo + 44] = morph_target_c;
}


// unpacks the tightly packed primitive-indices buffer from unsigned chars back to uint_32 
// needed, as HLSL not nativly supports single byte-sized variables
uint SampleTriangleBufferAsCharArray(uint i, uint mesh_id)
{
    uint array_pos = i / 4;
    uint inner_pos = i % 4;
    uint ret = primitiveIndicesBuffers[mesh_id][array_pos];
    ret = ret >> (inner_pos * 8);
    return ret & 0x000000ff;
}




uint Hash(uint x)
{
    x ^= x >> 16;
    x *= 0x85ebca6b; // Large prime
    x ^= x >> 13;
    x *= 0xc2b2ae35; // Another large prime
    x ^= x >> 16;
    return x;
}


// Generates a pseudo random float by hasing the provided seed
// Resulting float value lies between -1 and 1
float Random(uint seed)
{
    return float(Hash(seed * 3)) / 4294967295.0; // 2^32-1
}


// Returns a cheap approximation for a fully saturated Color value for the provided HUE
float4 Rainbow(float factor)
{
    float3 col = float3(abs(factor * 6.0 - 3.0) - 1.0, 2.0 - abs(factor * 6.0 - 2.0), 2.0 - abs(factor * 6.0 - 4.0));
    return float4(clamp(col, float3(0.0, 0.0, 0.0), float3(1.0, 1.0, 1.0)), 1.0);
}


// Returns the expected LOD on a vertex level while using the same distance function as its Task-Shader equivalent
// bounding sphere parameter has to be passed in world space coordiantes
// tessellation is represented as negative return value
float getExpectedLoDLevel(float4 position)
{
    float cam_dist = distance(constants.CameraWorldPos, position.xyz);
    return log2(cam_dist / constants.LoD_Scale);
}


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
    float4 xTex = heightMapTexture.SampleLevel(heightMapSampler, uvX, 0);
    float4 yTex = heightMapTexture.SampleLevel(heightMapSampler, uvY, 0);
    float4 zTex = heightMapTexture.SampleLevel(heightMapSampler, uvZ, 0);

    // Weighted blend
    float4 result = xTex * blend.x + yTex * blend.y + zTex * blend.z;
    //result = clamp(result, float4(0, 0, 0, 1), float4(1, 1, 1, 1));
    //return float4(1, 1, 1, 1);
    return result;
}


// mesh shader entry point
[shader("mesh")]
[numthreads(GROUP_SIZE, 1, 1)]
[outputtopology("triangle")]
void main(in uint I : SV_GroupIndex,
          in uint gid : SV_GroupID,
          in payload S_Payload gs_Payload,
          out indices uint3 tris[MAX_MESHLET_PRIMITIVE_COUNT],
          out vertices PixelShaderInput verts[MAX_MESHLET_VERTEX_COUNT])
{
    // retrieve a unique task from the global payload buffer for each individual mesh shader workgroup
    S_PayloadEntry payload_task = payloadBuffer[gid + gs_Payload.global_payload_offset];

    S_SceneObject scene_object = objectsBuffer[payload_task.object_id];
    S_Meshlet meshlet = meshletBuffers[scene_object.mesh_id][payload_task.meshlet_id];
    
    // optional tessellation
    bool do_tessellation = ((constants.BoolConstants & TRESSELLATION_BIT_POS) && payload_task.tessellation_grade);
    if (do_tessellation)
    {
        // Reset the group shared counter variables from the first thread in the group
        if (I == 0)
        {
            gs_triangle_count = 0;
            gs_vertex_count = 0;
        }
        GroupMemoryBarrierWithGroupSync();
        
        // tessellate each individual triangle
        for (uint p = payload_task.tessellation_triangle_offset + I; p < payload_task.tessellation_triangle_count + payload_task.tessellation_triangle_offset; p += GROUP_SIZE)
        {
            uint vertex_index_a = meshlet.vertex_offset + SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 0, scene_object.mesh_id);
            uint vertex_index_b = meshlet.vertex_offset + SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 1, scene_object.mesh_id);
            uint vertex_index_c = meshlet.vertex_offset + SampleTriangleBufferAsCharArray(meshlet.triangle_offset + p * 3 + 2, scene_object.mesh_id);
            
        
            if (payload_task.tessellation_grade == 3)
                lvl3Tessellation(uint3(vertex_index_a, vertex_index_b, vertex_index_c), scene_object.mesh_id);
            else if (payload_task.tessellation_grade == 2)
                lvl2Tessellation(uint3(vertex_index_a, vertex_index_b, vertex_index_c), scene_object.mesh_id);
            else
                lvl1Tessellation(uint3(vertex_index_a, vertex_index_b, vertex_index_c), scene_object.mesh_id);
        }
    }
    else
    {
        // Reset the group shared counter variables from the first thread in the group
        if (I == 0)
        {
            gs_triangle_count = meshlet.triangle_count;
            gs_vertex_count = meshlet.vertex_count;
        }
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    SetMeshOutputCounts(gs_vertex_count, gs_triangle_count);
    
    if (do_tessellation)
    {
        for (uint v = I; v < gs_vertex_count; v += GROUP_SIZE)
        {
            S_Vertex vertex = gs_vertices[v];
            
            float4 original_world_pos = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
            int current_vertex_lod = gs_tessellation_levels[v] * -1;
            float expected_lod = max(getExpectedLoDLevel(original_world_pos), MAX_TESSELLATION_LEVEL * -1);
            float lerp_value = max(expected_lod - current_vertex_lod, 0);
            
            // itterativly determine the correct vertex indices according to the vertex specific LOD
            // vertex LOD is generally lower than that of the whole meshlet. Geo-Morphing smootly blends between discrete levels 
            if ((constants.BoolConstants & GEO_MORPHING_BIT_POS) && !(constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS))
            {
                int vertex_index = v;
                // for vertex lod >= 0 the vertex_index variable actually represents a vertex-index-index and points into the regular vertex-indices / morph buffer
                while (lerp_value > 1)
                {
                    if (current_vertex_lod <= 0)
                    {
      
                        vertex_index = gs_tessellation_morph_targets[vertex_index];
                        current_vertex_lod++;

                        lerp_value--;
                    }
                    else
                    {
                        vertex_index = morphIndicesBuffers[scene_object.mesh_id][vertex_index];
                        current_vertex_lod++;
                        lerp_value--;
                    }
                }
                
                S_Vertex morph_target_vertex;
                // set both vertex and morph target according to current_lod case
                if (current_vertex_lod < 0)
                {
                    vertex = gs_vertices[vertex_index];
                    morph_target_vertex = gs_vertices[gs_tessellation_morph_targets[vertex_index]];
                }
                else if (current_vertex_lod == 0)
                {
                    vertex = gs_vertices[vertex_index];
                    morph_target_vertex = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][gs_tessellation_morph_targets[vertex_index]]];
                }
                else
                {
                    vertex = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][vertex_index]];
                    morph_target_vertex = verticesBuffers[scene_object.mesh_id][vertexIndicesBuffers[scene_object.mesh_id][morphIndicesBuffers[scene_object.mesh_id][vertex_index]]];
                }
                
                
                vertex = linearVertexInterpolation(vertex, morph_target_vertex, lerp_value);
            }
            
            
            // transform vertices from object to world to clip space
            // optional vertex displacement according to heightmap
            // populate outpiut buffers
            float4 w_pos = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
            float4 normal = mul(float4(vertex.normal.xyz, 0), scene_object.object_matrix);
            verts[v].WPos = w_pos.xyz;
            if (constants.BoolConstants & TRI_PLANAR_TEXTURE_MAPPING_BIT_POS)
                w_pos = float4(w_pos.xyz + normal.xyz * TriplanarSample(w_pos.xyz, normal.xyz).x * constants.HeightMapDisplacementScale, 1);
            else
                w_pos = float4(w_pos.xyz + normal.xyz * heightMapTexture.SampleLevel(heightMapSampler, vertex.uv.xy, 0).x * constants.HeightMapDisplacementScale, 1);
            verts[v].Pos = mul(w_pos, constants.ViewProjMat);
            verts[v].UV = vertex.uv.xy;
            verts[v].WNormal = normal.xyz;
            
            
            // set vertex color according to debug shading selection
            float brightness = clamp(clamp(dot(normalize(float3(1, 1, 1)), normal.xyz), 0, 1) + clamp(dot(normalize(float3(-2, 1, -1)), normal.xyz), 0, 0.6), 0.05, 1);
            if (!(constants.BoolConstants & NORMAL_LIGHTING_BIT_POS))
                brightness = 1;
            if (constants.shadingSelection == DEFAULT_SHADING)
            {
                verts[v].Color = brightness;
            }
            else if (constants.shadingSelection == MESHLETS_SHADING)
            {
                verts[v].Color = Rainbow(Random(gid + gs_Payload.global_payload_offset)) * brightness;
            }
            else if (constants.shadingSelection == LOD_SHADING)
            {
                if ((constants.BoolConstants & GEO_MORPHING_BIT_POS) && !(constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS))
                    verts[v].Color = Rainbow(frac((DEBUG_COLOR_SPREAD + max(expected_lod, -MAX_TESSELLATION_LEVEL) / DEBUG_COLOR_SPREAD))) * brightness;
                else
                    verts[v].Color = Rainbow(frac((DEBUG_COLOR_SPREAD - payload_task.tessellation_grade) / DEBUG_COLOR_SPREAD)) * brightness;
            }
            else if (constants.shadingSelection == MESHLET_GROUP_SHADING)
            {
                verts[v].Color = Rainbow(Random(meshlet.group_id)) * brightness;
            }
            else if (constants.shadingSelection == TESSELLATION_LEVEL_SHADING)
            {
                verts[v].Color = float4(1, 1, 1, 1) * brightness;
            }
            else
            {
            // should not be shown. resulting color is probably going to be overwritten by pixel shader
                verts[v].Color = float4(1, 1, 0, 1);
            }
        }
        
        
        for (uint p = I; p < gs_triangle_count; p += GROUP_SIZE)
        {
            tris[p] = gs_triangles[p];
        }
        
        return;
    }
    
    
    // process meshlet vertices
    const uint vertexLoops = (MAX_MESHLET_VERTEX_COUNT + GROUP_SIZE - 1) / GROUP_SIZE;
    for (uint v_loop = 0; v_loop < vertexLoops; v_loop++)
    {
        uint v = I + v_loop * GROUP_SIZE;
        v = min(v, meshlet.vertex_count - 1);
        
        int vertexIndex = vertexIndicesBuffers[scene_object.mesh_id][meshlet.vertex_offset + v];
        S_Vertex vertex = verticesBuffers[scene_object.mesh_id][vertexIndex];
        
        float4 original_world_pos = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
        
        float expected_LoD = max(getExpectedLoDLevel(original_world_pos), 0);
        float lerp_value = expected_LoD - meshlet.discrete_level_of_detail;
        //lerp_value = constants.DebugFloatSliderValue;
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
        
        
        float4 w_pos = mul(float4(vertex.position.xyz, 1.0), scene_object.object_matrix);
        float4 normal = mul(float4(vertex.normal.xyz, 0), scene_object.object_matrix);
        verts[v].WPos = w_pos.xyz;
        if (constants.BoolConstants & TRI_PLANAR_TEXTURE_MAPPING_BIT_POS)
            w_pos = float4(w_pos.xyz + normal.xyz * TriplanarSample(w_pos.xyz, normal.xyz).x * constants.HeightMapDisplacementScale, 1);
        else
            w_pos = float4(w_pos.xyz + normal.xyz * heightMapTexture.SampleLevel(heightMapSampler, vertex.uv.xy, 0).x * constants.HeightMapDisplacementScale, 1);
        verts[v].Pos = mul(w_pos, constants.ViewProjMat);
        verts[v].UV = vertex.uv.xy;
        verts[v].WNormal = normal.xyz;
        
        
        float brightness = clamp(clamp(dot(normalize(float3(1, 1, 1)), normal.xyz), 0, 1) + clamp(dot(normalize(float3(-2, 1, -1)), normal.xyz), 0, 0.6), 0.05, 1);
        if (!(constants.BoolConstants & NORMAL_LIGHTING_BIT_POS))
            brightness = 1;
        if (constants.shadingSelection == DEFAULT_SHADING)
        {
            verts[v].Color = brightness;
        }
        else if (constants.shadingSelection == MESHLETS_SHADING)
        {
            verts[v].Color = Rainbow(Random(gid + gs_Payload.global_payload_offset)) * brightness;
        }
        else if (constants.shadingSelection == LOD_SHADING)
        {
            if ((constants.BoolConstants & GEO_MORPHING_BIT_POS) && !(constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS))
                verts[v].Color = Rainbow(expected_LoD / 6.0) * brightness;
            else
                verts[v].Color = Rainbow(meshlet.discrete_level_of_detail / 6.0) * brightness;
            
        }
        else if (constants.shadingSelection == MESHLET_GROUP_SHADING)
        {
            verts[v].Color = Rainbow(Random(meshlet.group_id)) * brightness;
        }
        else if (constants.shadingSelection == TESSELLATION_LEVEL_SHADING)
        {
            verts[v].Color = float4(0.2, 0.2, 0.2, 1) * brightness;
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