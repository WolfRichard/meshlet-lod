#include "ViewDependentStructures.fxh"
#pragma optimize("", off)
//#pragma optimize("", on)


ConstantBuffer<S_Constants> constants                       : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer               : register(t0, space0);

RWStructuredBuffer<S_WorkQueueEntry> workQueue              : register(u0, space0);
RWStructuredBuffer<S_WorkQueueCounters> workQueueCounters   : register(u1, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]                : register(t0, space1);


// Payload will be used in the mesh shader.
groupshared S_Payload gs_Payload;

// The number of meshlets that are dispatched by this thread group,
groupshared uint gs_MeshletCount;




void appendTask(S_WorkQueueEntry new_task)
{
    uint task_index = 0;
    InterlockedAdd(workQueueCounters[0].tail, 1, task_index);    // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task;     // loop index because of ring buffer structure
}


bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint task_index, prev_index;
    do
    {
        prev_index = workQueueCounters[0].head; 
        uint end   = workQueueCounters[0].tail; 

        if ((prev_index >= end) || (gs_MeshletCount >= MAX_EMITTED_MESHLETS_PER_WORK_GROUP - GROUP_SIZE))
            return false; // No work available
        
        // Try to atomically update 'beginCounter' only if it hasn’t changed
        InterlockedCompareExchange(workQueueCounters[0].head, prev_index, prev_index + 1, task_index);
        
    } while ((task_index != prev_index)); // lost race condition in between checking for available work and incrementing counter --> so check if there is other work available to grab
    
    // Load the task safely as thread won race condition after checking that there was still work to do
    out_task = workQueue[task_index % WORK_QUEUE_SIZE];
    return true;
}




// Frustum culling in world space
bool isInFrustum(float3 center, float radius)
{
    // TODO: remove early termination 
    return true;
    // TODO: FIX CULLING currently cuts away to early
    float4 f4Center = float4(center, 1.0);
    for (int i = 0; i < 6; ++i)
    {
        if (dot(constants.Frustum[i], f4Center) < -radius)
            return false;
    }
    return true;
}


float ExtractMaxScaleFactor(float4x4 m)
{
    float sx = length(m[0].xyz);
    float sy = length(m[1].xyz);
    float sz = length(m[2].xyz);
   
    return max(sx, max(sy, sz));
}

float getExpectedLoDLevel(S_BoundingSphere bounding_sphere)// bounding sphere must be in world space!
{
    float cam_dist = max(distance(constants.CameraWorldPos, bounding_sphere.center) - bounding_sphere.radius, 0);
    //float cam_dist = distance(constants.CameraWorldPos, bounding_sphere.center);
    return max(log2(cam_dist / constants.LoD_Scale), 0);
}


bool groupSimplificationIsPreciseEnough(S_BoundingSphere bounding_sphere, uint lod_level) // bounding sphere must be in world space!
{
    //return lod_level <= 0;
    return lod_level <= getExpectedLoDLevel(bounding_sphere);
}

float getScreenSpaceErrorInPixels(S_BoundingSphere bounding_sphere) // bounding sphere must be in clip space
{
    float d2 = dot(bounding_sphere.center, bounding_sphere.center);
    float r2 = bounding_sphere.radius * bounding_sphere.radius;
    float sphere_diameter_uv = max(constants.ProjMat[0][0], constants.ProjMat[1][1]) * bounding_sphere.radius / sqrt(d2 - r2);
    float view_size = max(constants.ScreenWidth, constants.ScreenHeight);
    return sphere_diameter_uv * view_size;
}

bool isPreciseEnough(S_BoundingSphere bounding_sphere) // bounding sphere must be in clip space
{
    return getScreenSpaceErrorInPixels(bounding_sphere) < constants.LoD_Scale; //    1.0;
}

void delayBeforeNextWorkQueueCheck(uint number_of_attemps_scince_last_work)
{
    float dummy_work_load;
    for (uint i = 0; i < number_of_attemps_scince_last_work * 1; i++)
        dummy_work_load = sqrt(max(i, 2));
}

void queueMeshletForDispatch(uint meshlet_index, uint object_index)
{
    S_SceneObject scene_object = objectsBuffer[object_index];
    S_BoundingSphere world_space_bounding_sphere;
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][meshlet_index];
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz;
    world_space_bounding_sphere.radius = current_meshlet.bounding_sphere.radius * ExtractMaxScaleFactor(scene_object.object_matrix);
    
        uint payload_index = 0;
        InterlockedAdd(gs_MeshletCount, 1, payload_index); // Undefined behavior if max meshlet count per work group is set to low 
        gs_Payload.tasks[payload_index].meshlet_id = meshlet_index;
        gs_Payload.tasks[payload_index].object_id = object_index;
    
}

void processTask(S_WorkQueueEntry task)
{
    S_SceneObject scene_object = objectsBuffer[task.scene_object_id];
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][task.meshlet_id];
        
    bool precise_enough;
            
    S_BoundingSphere world_space_bounding_sphere;
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz; // object --> world space
    world_space_bounding_sphere.radius = current_meshlet.base_error * ExtractMaxScaleFactor(scene_object.object_matrix);
        
    if (constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS)
    {
        world_space_bounding_sphere.center = mul(float4(world_space_bounding_sphere.center, 1), constants.ViewMat).xyz; // world --> clip space
        precise_enough = isPreciseEnough(world_space_bounding_sphere); // actually in screen clip space
    }
    else
    {
        precise_enough = groupSimplificationIsPreciseEnough(world_space_bounding_sphere, current_meshlet.discrete_level_of_detail);
    }
        
    if (precise_enough)
    {
        if (isInFrustum(world_space_bounding_sphere.center, world_space_bounding_sphere.radius))
            queueMeshletForDispatch(task.meshlet_id, task.scene_object_id);
    }
    else // append all child meshlets into work queue if current simplification isnt precise enough (zero children when current meshlet is a leaf node)
    {
        for (uint child_meshlet_index_index = 0; child_meshlet_index_index < current_meshlet.child_count; child_meshlet_index_index++)
        {
            S_WorkQueueEntry new_task;
            new_task.scene_object_id = task.scene_object_id;
            new_task.meshlet_id = current_meshlet.child_meshlets[child_meshlet_index_index];
            appendTask(new_task);
        }
    }
}


[numthreads(GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
    uint global_thread_index = I + wg * GROUP_SIZE; 
    
   // Reset the group shared counter variables from the first thread in the group
    if (I == 0)
    {
        gs_MeshletCount = 0;
    }
    GroupMemoryBarrierWithGroupSync();
 
    // object culling and extracting of root nodes into work queue
    for (uint scene_object_index = global_thread_index; scene_object_index < constants.SceneObjectCount; scene_object_index += PERSISTENT_THREAD_COUNT)
    {
        S_SceneObject scene_object = objectsBuffer[scene_object_index];
        
        if (isInFrustum(scene_object.bounding_sphere.center, scene_object.bounding_sphere.radius))
        {
            for (uint i = 0; i < GROUP_SPLIT_COUNT; i++)
            {
                S_WorkQueueEntry new_task;
                new_task.meshlet_id = meshletBuffers[scene_object.mesh_id][0].child_meshlets[i];
                new_task.scene_object_id = scene_object_index;
                appendTask(new_task);
            }
        }
    }
        
    AllMemoryBarrierWithGroupSync();
       
    // processing of work queue entries    
    S_WorkQueueEntry current_task;
    // check if new task have been added aslong as there are still tasks beeing processed by other threads
    //while (consumeTask(current_task))
    bool keepSpinLock = ((workQueueCounters[0].done_processing < workQueueCounters[0].tail) && (gs_MeshletCount < MAX_EMITTED_MESHLETS_PER_WORK_GROUP - GROUP_SIZE));
    
    while (keepSpinLock)
    {
        if (consumeTask(current_task))
        {
            processTask(current_task);
            uint last_task_index = 0;
            InterlockedAdd(workQueueCounters[0].done_processing, 1, last_task_index);
        }
            
        keepSpinLock = ((workQueueCounters[0].done_processing < workQueueCounters[0].tail) && (gs_MeshletCount < MAX_EMITTED_MESHLETS_PER_WORK_GROUP - GROUP_SIZE));
    }     
    //after work queue is empty dispatch all collected meshlets
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
}
