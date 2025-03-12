# pragma once


#define GROUP_SIZE 32

//#define MAX_MESHLET_VERTEX_COUNT 64
#define MAX_MESHLET_VERTEX_COUNT 128

//#define MAX_MESHLET_PRIMITIVE_COUNT 124
#define MAX_MESHLET_PRIMITIVE_COUNT 254

#define PI 3.1415927

struct S_CullingInfo // culling data is presented in object space
{
    float3 bounding_sphere_center;
    float bounding_sphere_radius;

    float3 normal_cone_apex;
    float normal_cone_cutoff;
    float3 normal_cone_axis;

    float byte_alignement; // has no use other than ensuring that the structure is 16-byte aligned
};

struct S_DrawTask
{
    uint vertex_count;
    uint triangle_count;
    uint vertex_offset;
    uint triangle_offset;

    S_CullingInfo culling_info;
};

// bit positions for constant bools
#define FRUSTUM_CULLING_BIT_POS          1 // 00000000000000000000000000000001
#define CONE_CULLING_BIT_POS             2 // 00000000000000000000000000000010
#define ENABLE_DEBUG_VISUALS_BIT_POS     4 // 00000000000000000000000000000100
#define DEBUG_MESHLETS                  16 // 00000000000000000000000000001000
#define DEBUG_BONES                     32 // 00000000000000000000000000010000
#define DEBUG_LOD                       64 // 00000000000000000000000000100000
#define ENABLE_OBJECT_LOD              128 // 00000000000000000000000001000000
#define ENABLE_MESHLET_LOD             256 // 00000000000000000000000010000000

struct S_Constants
{
    float4x4 ViewProjMat;
    float4x4 ViewMat;
    float4 Frustum[6];
    
    float3 CameraWorldPos;
    float CoTanHalfFoV;
    
    float CurrTime;
    
    uint BoolConstants; // can hold up to 32 booleans but only uses the space of a single word in the root signature
};


// Payload from task shader stage to mesh shader
// Payload size must be less than 16kb.
struct S_Payload
{
    uint vertex_count[GROUP_SIZE];
    uint triangle_count[GROUP_SIZE];
    
    uint vertex_offset[GROUP_SIZE];
    uint triangle_offset[GROUP_SIZE];
    
    uint meshlet_LoD[GROUP_SIZE];
    
    uint object_id;
    
};

struct S_Vertex
{
    float4 position;
    float4 uv;
    float4 normal;
    float4 color;

};

struct S_SceneObject
{
    float4x4 object_matrix;
    
    float3 bounding_sphere_center; // in world space
    float bounding_sphere_radius; // in world space
    
    uint mesh_id;
    
    float3 byte_allignement; // used to keep 16 byte allignement for structured buffer
};


struct S_CommandStructure
{
    uint instanceID;
    float level_of_detail;
    uint3 dispatchArguments;
};