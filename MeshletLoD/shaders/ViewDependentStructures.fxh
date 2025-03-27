# pragma once


#define GROUP_SIZE 32

#define PERSISTENT_THREAD_COUNT 1024

#define WORK_QUEUE_SIZE 65536           // (2^16) Workqueue is implemented as ring-buffer, 
                                        //queue size allows correct indexing and should be able to contain the maximum of simultanious queue tasks

#define MAX_EMITTED_MESHLETS_PER_WORK_GROUP 1024 // Maximum of 32 Meshlets per Task Shader Thread (limited by maximum size of payload)

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
    S_BoundingSphere bounding_sphere; // in object space
    uint vertex_count;      
    uint triangle_count;    // how many primitives does the meshlet consists of
    uint vertex_offset;
    uint triangle_offset;   // offset to the first triangle of the meshlet (offset is counted in indices not triangles !!!)
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
#define FRUSTUM_CULLING_BIT_POS          1  // 00000000000000000000000000000001
//#define CONE_CULLING_BIT_POS             2  // 00000000000000000000000000000010
//#define ENABLE_DEBUG_VISUALS_BIT_POS     4  // 00000000000000000000000000000100
//#define DEBUG_MESHLETS                  16  // 00000000000000000000000000001000
//#define DEBUG_BONES                     32  // 00000000000000000000000000010000
//#define DEBUG_LOD                       64  // 00000000000000000000000000100000
//#define ENABLE_OBJECT_LOD              128  // 00000000000000000000000001000000
//#define ENABLE_MESHLET_LOD             256  // 00000000000000000000000010000000


enum ShadingMode
{
    DEFAULT_SHADING,
    DEBUG_MESHLET_SHADING,
    DEBUG_LOD_SHADING,
    DEBUG_WORLD_POS,
    DEBUG_MESHLET_GROUP
};


struct S_Constants
{
    float4x4 ViewProjMat;
    float4x4 ViewMat;
    float4 Frustum[6];
    
    float3 CameraWorldPos;
    float CoTanHalfFoV;
    float LoD_Scale;
    
    float CurrTime;
    
    float MaxScaleFactor_ViewProjMat;
    
    uint SceneObjectCount;
    
    ShadingMode shadingSelection;
    uint BoolConstants;                     // can hold up to 32 booleans but only uses the space of a single word in the root signature
};


// every mesh shader thread group processes one payload entry which translates to one meshlet 
struct S_PayloadEntry
{
    uint meshlet_id;
    uint object_id;
    float lod_morphing; // 0 = vertices stay the same; 1 = vertices are completely pushed towards their parent vertex
    uint lod_tree_depth;
};

// Payload from task shader stage to mesh shader
// Payload size must be less than 16kb.
struct S_Payload
{
    S_PayloadEntry tasks[MAX_EMITTED_MESHLETS_PER_WORK_GROUP];
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
    
    uint root_group_id; // could also be stored in seperate "mesh" buffer to avoid storing this multiple times per instance
    uint mesh_group_count; // could also be stored in seperate "mesh" buffer to avoid storing this multiple times per instance
    uint mesh_meshlet_count; // could also be stored in seperate "mesh" buffer to avoid storing this multiple times per instance
    
    
    float2 byte_allignement;            // used to keep 16 byte allignement for structured buffer
};

struct S_WorkQueueEntry
{
    uint scene_object_id;
    uint group_id;
};
