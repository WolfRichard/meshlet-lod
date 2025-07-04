# pragma once

#define GROUP_SIZE 32 // has to stay at 32 to match wave size, also not allowed to be multiple of 32 as tis would break the task shader wave intrinsics logic

#define PERSISTENT_THREAD_COUNT 81920    // has to be multiple of GROUP_SIZE

#define WORK_QUEUE_SIZE 4194304         // (2^22) Workqueue is implemented as ring-buffer, queue size allows correct indexing and should be able to contain the maximum of simultanious queue tasks

#define MAX_DISPATCH_MESH_GROUP_COUNT 4194304 // 2^22 https://microsoft.github.io/DirectX-Specs/d3d/MeshShader.html#dispatchmesh-api

// tessellation constants
#define MAX_TESSELLATION_LEVEL 3
// look up table for triangle and vertex counts per tessellation level (from lvl 1-3)
static const uint TESSELLATION_INITIAL_TRIANGLE_COUNTS[3]   = {  40,  15,   3 };
//static const uint TESSELLATION_RESULTING_TRIANGLE_COUNTS[3] = { 120, 240, 192 };
//static const uint TESSELLATION_INITIAL_VERTEX_COUNTS[3]     = { 120,  45,   9 };
//static const uint TESSELLATION_RESULTING_VERTEX_COUNTS[3]   = { 240, 225, 135 };

#define DEBUG_COLOR_SPREAD 5.0 // the number of individual colors used for debugging 

#define MAX_MESHLET_VERTEX_COUNT 254

#define MAX_MESHLET_PRIMITIVE_COUNT 254

#define PI 3.1415927
#ifndef FLT_MAX
    #define FLT_MAX 3.402823466e+38
#endif

struct S_BoundingSphere
{
    float3 center;
    float radius;
};

#define GROUP_MERGE_COUNT 4 // number of meshlets that are grouped together
#define GROUP_SPLIT_COUNT 2 // number of meshlets that are generated out of the merged meshlets

struct S_Meshlet
{
    S_BoundingSphere bounding_sphere; // in object space
    uint vertex_count;      
    uint triangle_count;              // how many primitives does the meshlet consists of
    uint vertex_offset;
    uint triangle_offset;             // offset to the first triangle of the meshlet (offset is counted in indices not triangles !!!)
    
    S_BoundingSphere simplified_group_bounds;
    float base_error;                       // if the meshlet is a sub-mesh of the original geometry this is set to 0
    float simplification_error;             // (parent error) if meshlet has no parent because its part of the root then this is set to FLOAT_MAX
    uint discrete_level_of_detail;          // the discrete level of innacuracy for the meshlet (higher number means less precise mesh resolution)
    uint group_id;                          // currently just for debug-rendering
    
    uint child_meshlets[GROUP_MERGE_COUNT]; // indices that point into the meshlet buffer towards the child meshlets that were used during group simplification to produce the simplififed geometry of the current meshlet
    
    uint child_count;
    float3 byte_allignement;                // used to keep 16 byte allignement for structured buffer
};


// node of the binary tree, each group consists of 2 Meshlets
struct S_MeshletGroup
{
    S_BoundingSphere bounding_sphere;               // in object space, used for LoD-decision, radius = simpification error so maxiumum Offset that a vertex was moved from its original position during simplification in obect space. (equals 0 for leaf nodes, should equal FLOAT_MAX for parent of root node)
    int parent;                                     // index to parent node, -1 if root
    int children[GROUP_SPLIT_COUNT * 2];            // indices to child nodes -1 if no child (* 2 as buffer for partion imperfections )
    uint childCount;                                // can also be used as bool to check if a node is a leaf
    uint simplified_meshlets[GROUP_SPLIT_COUNT];    // indices to the simplified meshlets that are part of this group 
    uint meshlets[GROUP_MERGE_COUNT * 2];           // indices to the more granular meshlets that are part of this group (0 is a dummy meshlet with no content) (* 2 as buffer for partion imperfections)       
    uint meshlet_count;                             // allows up to merge count * 2 meshlets, this is necessary because graph partitioning wont always create perfect groups of 4
    uint hierarchy_tree_depth;
    
    float2 byte_allignement;                        // used to keep 16 byte allignement for structured buffer
};


// bit positions for constant bools
#define FRUSTUM_CULLING_BIT_POS                 1  // 00000000000000000000000000000001
#define SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS    2  // 00000000000000000000000000000010
#define GEO_MORPHING_BIT_POS                    4  // 00000000000000000000000000000100
#define TRESSELLATION_BIT_POS                  16  // 00000000000000000000000000001000
#define TRI_PLANAR_TEXTURE_MAPPING_BIT_POS     32  // 00000000000000000000000000010000
#define NORMAL_LIGHTING_BIT_POS                64  // 00000000000000000000000000100000
//#define PLACEHOLDER_BIT_POS                 128  // 00000000000000000000000001000000
//#define PLACEHOLDER_BIT_POS                 256  // 00000000000000000000000010000000
//#define PLACEHOLDER_BIT_POS                 512  // 00000000000000000000000100000000
//#define PLACEHOLDER_BIT_POS                1024  // 00000000000000000000001000000000
//#define PLACEHOLDER_BIT_POS                2048  // 00000000000000000000010000000000
//#define PLACEHOLDER_BIT_POS                4096  // 00000000000000000000100000000000
//#define PLACEHOLDER_BIT_POS                8192  // 00000000000000000001000000000000


enum ShadingMode
{
    DEFAULT_SHADING,
    MESHLETS_SHADING,
    MESHLET_GROUP_SHADING,
    LOD_SHADING,
    TESSELLATION_LEVEL_SHADING,
    HEIGHT_MAP_SHADING,
};


// defines constant buffer structure
struct S_Constants
{
    float4x4 ViewProjMat;
    float4x4 ViewMat;
    float4x4 ProjMat;
    float4 Frustum[6];
    
    float3 CameraWorldPos;
    float CoTanHalfFoV;
    float LoD_Scale;
    
    float CurrTime;
    uint ScreenWidth; // in pixels
    uint ScreenHeight; // in pixels
    
    float MaxScaleFactor_ViewProjMat;
    
    float DebugFloatSliderValue;
    
    float TriPlanarMappingScale;
    float TriPlanarBlendGrade;
    
    float HeightMapDisplacementScale;
    
    uint SceneObjectCount;
    
    ShadingMode shadingSelection;
    uint BoolConstants;                     // can hold up to 32 booleans but only uses the space of a single word in the root signature
};


// every mesh shader thread group processes one payload entry which translates to one meshlet 
struct S_PayloadEntry
{
    uint meshlet_id;
    uint object_id;
    uint tessellation_grade;
    uint tessellation_triangle_offset;
    uint tessellation_triangle_count;
};


// Payload from task shader stage to mesh shader
// Payload size must be less than 16kb.
struct S_Payload
{
    uint global_payload_offset;
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
    
    uint mesh_meshlet_count;            // could also be stored in seperate "mesh" buffer to avoid storing this multiple times per instance
    
    float2 byte_allignement;            // used to keep 16 byte allignement for structured buffer
};

struct S_WorkQueueEntry
{
    uint scene_object_id;
    uint meshlet_id;
};

struct S_WorkQueueCounters
{
    uint work_queue_head; // index to first element of the FIFO ring buffer work queue (if head == tail then queue is empty)
    uint work_queue_tail; // index of the next free available space where the next element could be added (tail is not part of the queue)
    uint payload_head;    // index from what point the meshlets have not been emitted yet
    uint payload_tail;    // index to the next available free meshlet slot
};

