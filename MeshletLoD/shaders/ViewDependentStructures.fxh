# pragma once


#define GROUP_SIZE 32

//#define MAX_MESHLET_VERTEX_COUNT 64
#define MAX_MESHLET_VERTEX_COUNT 128

//#define MAX_MESHLET_PRIMITIVE_COUNT 124
#define MAX_MESHLET_PRIMITIVE_COUNT 254

#define PI 3.1415927

struct S_BoundingSphere
{
    float3 center;
    float radius;
};

#define GROUP_MERGE_COUNT 4 // number of meshlets that are grouped together
#define GROUP_SPLIT_COUNT 2 // number of meshlets that are generated out of the merged meshlets

struct S_Meshlet
{
    uint vertex_count;      
    uint triangle_count;    // how many primitives does the meshlet consists of
    uint vertex_offset;
    uint triangle_offset;   // offset to the first triangle of the meshlet (offset is counted in indices not triangles !!!)
};


// node of the binary tree, each group consists of 2 Meshlets
struct S_Meshlet_Group
{
    S_BoundingSphere bounding_sphere;           // in object space
    int parent;                                 // index to parent node, -1 if root
    int children[2];                            // indices to child nodes -1 if no child
    int isLeaf;                                 // = 0 means has children, != 0 means no children
    int simplified_meshlets[GROUP_SPLIT_COUNT]; // indices to the simplified meshlets that are part of this group
    int meshlets[GROUP_MERGE_COUNT];            // indices to the more granular meshlets that are part of this group 
    
    float2 byte_allignement;            // used to keep 16 byte allignement for structured buffer
};


// bit positions for constant bools
#define FRUSTUM_CULLING_BIT_POS          1  // 00000000000000000000000000000001
#define CONE_CULLING_BIT_POS             2  // 00000000000000000000000000000010
#define ENABLE_DEBUG_VISUALS_BIT_POS     4  // 00000000000000000000000000000100
#define DEBUG_MESHLETS                  16  // 00000000000000000000000000001000
#define DEBUG_BONES                     32  // 00000000000000000000000000010000
#define DEBUG_LOD                       64  // 00000000000000000000000000100000
#define ENABLE_OBJECT_LOD              128  // 00000000000000000000000001000000
#define ENABLE_MESHLET_LOD             256  // 00000000000000000000000010000000

struct S_Constants
{
    float4x4 ViewProjMat;
    float4x4 ViewMat;
    float4 Frustum[6];
    
    float3 CameraWorldPos;
    float CoTanHalfFoV;
    
    float CurrTime;
    
    uint BoolConstants;                     // can hold up to 32 booleans but only uses the space of a single word in the root signature
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
    float4x4 object_matrix;             // transforms from objcet to world space
    
    S_BoundingSphere bounding_sphere;   // in world space
    
    uint mesh_id;                       // the mesh of the object
    
    float3 byte_allignement;            // used to keep 16 byte allignement for structured buffer
};

