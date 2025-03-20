#include "ViewDependentStructures.fxh"



StructuredBuffer<S_SceneObject> objectsBuffer       : register(t0, space0);
ConstantBuffer<S_Constants> constants               : register(b0, space0);

RWStructuredBuffer<S_WorkQueueEntry> workQueue      : register(u0, space0);

RWStructuredBuffer<S_PayloadEntry> globalPayload    : register(u1, space0);


// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]        : register(t0, space1);
StructuredBuffer<S_MeshletGroup> groupBuffers[]          : register(t0, space2);


// Workqueue read- & write-Indices
groupshared uint gs_WorkQueueHead;
groupshared uint gs_WorkQueueTail;

void appendTask(S_WorkQueueEntry new_task)
{
    uint task_index;
    InterlockedAdd(gs_WorkQueueTail, 1, task_index); // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task; // loop index for ring buffer structure
}

bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint task_index;
    InterlockedAdd(gs_WorkQueueHead, 1, task_index); // Atomically increment beginCounter

    // Check if the queue is empty
    if (task_index >= gs_WorkQueueTail) // Compare with endCounter
        return false;

    out_task = workQueue[task_index % WORK_QUEUE_SIZE]; // Circular buffer
    return true;
}


// Payload will be used in the mesh shader.
groupshared S_Payload gs_Payload;

// The number of meshlets that are visible by the camera,
groupshared uint gs_MeshletCount;





// Frustum Culling
bool isInFrustum(float3 center, float radius)
{
    float4 f4Center = float4(center, 1.0);
    for (int i = 0; i < 6; ++i)
    {
        if (dot(constants.Frustum[i], f4Center) < -radius)
            return false;
    }
    return true;
}

// Cone Culling
bool IsShowingBackside(float3 coneApex, float3 coneAxis, float coneCutoff)
{
    float3 directionConeCamera = normalize(constants.CameraWorldPos - coneApex);
    return (dot(directionConeCamera, coneAxis) >= coneCutoff);
}

float CalcDetailLevel(float3 center, float radius)
{
    float3 pos = mul(float4(center, 1.0), constants.ViewMat).xyz;
    float dist2 = dot(pos, pos);
    
    // calculate the sphere size in screen space
    float size = 2 * constants.CoTanHalfFoV * radius / sqrt(dist2 - radius * radius);
    

    float level = clamp(1.0 - size, 0.0, 0.999);
    return level * level;
}

[numthreads(MAXIMUM_GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
   // Reset the counter variables from the first thread in the group
    if (I == 0)
    {
        gs_MeshletCount = 0;
        gs_WorkQueueHead = 0;
        gs_WorkQueueTail = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    // object culling and extracting of root nodes into work queue
    for (uint scene_object_index = I; scene_object_index < constants.SceneObjectCount; scene_object_index += MAXIMUM_GROUP_SIZE)
    {
        S_SceneObject scene_object = objectsBuffer[scene_object_index];
        
        if (isInFrustum(scene_object.bounding_sphere.center, scene_object.bounding_sphere.radius))
        {
            S_WorkQueueEntry new_task;
            new_task.group_id = scene_object.root_group_id;
            new_task.scene_object_id = scene_object_index;
        }
    }
    
    // processing of work queue entries
    S_WorkQueueEntry current_task;
    while (consumeTask(current_task))
    {
        
        S_SceneObject scene_object = objectsBuffer[current_task.scene_object_id];
        S_MeshletGroup current_group = groupBuffers[scene_object.mesh_id][current_task.group_id];
        
        
        
        // TODO: Meshlet Culling
        // TODO: LoD Decision
        // TODO: Tree Traversal
        // TODO: WorkQueue Append Group Children
        // TODO: increment gs_MeshletCount
    }
    
    
    // dispatch of mesh shader instances to process the
    if (I == 0)
    {
        gs_Payload.meshlet_count = gs_MeshletCount;
    }
    GroupMemoryBarrierWithGroupSync();
    
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
}
