# pragma once

#ifndef GROUP_SIZE
#    define GROUP_SIZE 32
#endif

#ifndef MAX_MESHLET_VERTEX_COUNT
#    define MAX_MESHLET_VERTEX_COUNT 64
#endif

#ifndef MAX_MESHLET_PRIMITIVE_COUNT
#    define MAX_MESHLET_PRIMITIVE_COUNT 124
#endif

#define PI 3.1415927

struct CullingInfo
{
    float3 bounding_sphere_center;
    float bounding_sphere_radius;

    float3 normal_cone_apex;
    float normal_cone_cutoff;
    float3 normal_cone_axis;

    float byte_alignement; // has no use other than ensuring that the structure is 16-byte aligned
};

struct DrawTask
{
    uint vertex_count;
    uint triangle_count;

    uint vertex_offset;
    uint triangle_offset;

    CullingInfo culling_info;
};

// bit positions for constant bools
#define FRUSTUM_CULLING_BIT_POS 1 // 00000000000000000000000000000001
#define CONE_CULLING_BIT_POS    2 // 00000000000000000000000000000010
#define DEBUG_VISUALS_BIT_POS   4 // 00000000000000000000000000000100

struct Constants
{
    float4x4 ViewProjMat;
    float4 Frustum[6];
    
    float3 CameraWorldPos;

    float CurrTime;
    
    uint BoolConstants; // can hold up to 32 booleans but only uses the space of a single word in the root signature
};


// Payload from task shader stage to mesh shader
// Payload size must be less than 16kb.
struct Payload
{
    uint vertex_count[GROUP_SIZE];
    uint triangle_count[GROUP_SIZE];
    
    uint vertex_offset[GROUP_SIZE];
    uint triangle_offset[GROUP_SIZE];
    
    float3 allignement[GROUP_SIZE];
    
    uint object_id[GROUP_SIZE];
};

struct CustomVertex
{
    float4 position;
    float4 uv;
    float4 normal;
    float4 color;
};

struct IndirectMeshDrawArguments
{
    uint task_count;
    uint first_task;
};

struct SceneObject
{
    float4x4 object_matrix;
    
    float3 bounding_sphere_center;
    float bounding_sphere_radius;
    
    uint mesh_id;
    
    float3 byte_allignement; // used to keep 16 byte allignement for structured buffer
};