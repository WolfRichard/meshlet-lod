#include "structures.fxh"


// Objects Buffer
StructuredBuffer<SceneObject> objectsBuffer     : register(t0, space0);
// Per Mesh: Meshlet Counts Buffer
StructuredBuffer<uint> meshletCountsBuffer      : register(t1, space0);
// Draw task arguments
StructuredBuffer<DrawTask> drawTaskBuffers[]    : register(t0, space4);


cbuffer InstanceIDBuffer                        : register(b0, space0)
{
    Constants constantsBuffer;
};


cbuffer InstanceIDBuffer                        : register(b1, space0)
{
    uint object_id;
};

// Payload will be used in the mesh shader.
groupshared Payload s_Payload;

// The number of meshlets that are visible by the camera,
groupshared uint s_TaskCount;


// Frustum Culling
bool IsInFrustum(float3 center, float radius)
{
    float4 f4Center = float4(center, 1.0);
    for (int i = 0; i < 6; ++i)
    {
        if (dot(constantsBuffer.Frustum[i], f4Center) < -radius)
            return false;
    }
    return true;
}

// Cone Culling
bool IsShowingBackside(float3 coneApex, float3 coneAxis, float coneCutoff)
{
    return (dot(normalize(constantsBuffer.CameraWorldPos - coneApex), coneAxis) >= coneCutoff); // cos(acos(coneCutoff) * 2));

}


[numthreads(GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
    // Reset the counter from the first thread in the group
    if (I == 0)
    {
        s_TaskCount = 0;
    }

    GroupMemoryBarrierWithGroupSync();

    const uint gid = wg * GROUP_SIZE + I;
    
    
    
    //uint object_ID = 0; // incrementingConstant;
    
    if (gid < meshletCountsBuffer[objectsBuffer[object_id].mesh_id])
    {
        
        DrawTask task = drawTaskBuffers[objectsBuffer[object_id].mesh_id][gid];
        
        float3 bounding_sphere_center = mul(float4(task.culling_info.bounding_sphere_center, 1.0), objectsBuffer[object_id].object_matrix).xyz;
        float bounding_sphere_radius = length(mul(float4(task.culling_info.bounding_sphere_center + float3(task.culling_info.bounding_sphere_radius, 0, 0), 1.0), objectsBuffer[object_id].object_matrix).xyz - bounding_sphere_center);
        
        
        if ((!(constantsBuffer.BoolConstants & FRUSTUM_CULLING_BIT_POS) || IsInFrustum(bounding_sphere_center, bounding_sphere_radius)) &&
            (!(constantsBuffer.BoolConstants & CONE_CULLING_BIT_POS) || !IsShowingBackside(task.culling_info.normal_cone_apex, task.culling_info.normal_cone_axis, task.culling_info.normal_cone_cutoff)))
        {
        // Acquire an index that will be used to safely access the payload.
        // Each thread gets a unique index.
            uint index = 0;
            InterlockedAdd(s_TaskCount, 1, index);
    
            s_Payload.vertex_count[index] = task.vertex_count;
            s_Payload.triangle_count[index] = task.triangle_count;
    
            s_Payload.vertex_offset[index] = task.vertex_offset;
            s_Payload.triangle_offset[index] = task.triangle_offset;
            
            s_Payload.object_id[index] = object_id;
            //s_Payload.allignement[index] = float3(0, 0, 0);
        }
    
        
    }
   
    
   
    // All threads must complete their work so that we can read s_TaskCount
    GroupMemoryBarrierWithGroupSync();

    DispatchMesh(s_TaskCount, 1, 1, s_Payload);
    
    
    
    
    
    
    
    
    
}