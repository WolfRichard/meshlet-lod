#include "ViewDependentStructures.fxh"


ConstantBuffer<S_Constants> constants               : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer       : register(t0, space0);

RWStructuredBuffer<S_WorkQueueEntry> workQueue      : register(u0, space0);
// Workqueue read- & write-Indices
// 0  -> Begin Counter (read index)
// 1  -> End Counter (write index)
RWStructuredBuffer<uint> workQueueCounters          : register(u1, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]        : register(t0, space1);
StructuredBuffer<S_MeshletGroup> groupBuffers[]     : register(t0, space5);




void appendTask(S_WorkQueueEntry new_task)
{
    uint task_index;
    InterlockedAdd(workQueueCounters[1], 1, task_index);    // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task;     // loop index because of ring buffer structure
}

bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint task_index;
    InterlockedAdd(workQueueCounters[0], 1, task_index);    // Atomically increment beginCounter

    // Check if the queue is empty
    if (task_index >= workQueueCounters[1])                 // Compare with endCounter
        return false;

    out_task = workQueue[task_index % WORK_QUEUE_SIZE];     // loop index because of ring buffer structure
    return true;
}


// Payload will be used in the mesh shader.
groupshared S_Payload gs_Payload;

// The number of meshlets that are dispatched by this thread group,
groupshared uint gs_MeshletCount;


// Frustum culling in world space
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

bool groupSimplificationIsPreciseEnough(S_BoundingSphere bounding_sphere) // bounding sphere must be in world space!
{
    float3 pos = mul(float4(bounding_sphere.center, 1.0), constants.ViewMat).xyz;
    float dist2 = dot(pos, pos);
    
    // calculate the sphere size in screen space
    float size = constants.CoTanHalfFoV * bounding_sphere.radius / sqrt(dist2 - bounding_sphere.radius * bounding_sphere.radius);
    
    return (size < 1 / 480); // screen pixel count over FoV !!!!!! TODO: should be set in constants buffer !!!!!!!
}


void queueMeshletForDispatch(uint meshlet_index, uint object_index, float lod_blend_value)
{
    S_SceneObject scene_object = objectsBuffer[object_index];
    S_BoundingSphere world_space_bounding_sphere;
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][meshlet_index];
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz;
    world_space_bounding_sphere.radius = length(mul(float4(current_meshlet.bounding_sphere.radius, 0, 0, 0), scene_object.object_matrix).xyz); // theoretically only considers scale in X-Axes so technically incorrect
    if (isInFrustum(world_space_bounding_sphere.center, world_space_bounding_sphere.radius))
    {
        uint payload_index = 0;
        InterlockedAdd(gs_MeshletCount, 1, payload_index); // Undefined behavior if max meshlet count per work group is set to low 
        gs_Payload.tasks[payload_index].meshlet_id = meshlet_index;
        gs_Payload.tasks[payload_index].object_id = object_index;
        gs_Payload.tasks[payload_index].lod_morphing = lod_blend_value;
    }
}

[numthreads(GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
   // Reset the counter variables from the first thread in the group
    if (I == 0)
    {
        gs_MeshletCount = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    // object culling and extracting of root nodes into work queue
    for (uint scene_object_index = I; scene_object_index < constants.SceneObjectCount; scene_object_index += PERSISTENT_THREAD_COUNT)
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
        S_BoundingSphere world_space_bounding_sphere;
        world_space_bounding_sphere.center = mul(float4(current_group.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz;
        world_space_bounding_sphere.radius = length(mul(float4(current_group.bounding_sphere.radius, 0, 0, 0), scene_object.object_matrix).xyz); // theoretically only considers scale in X-Axes so technically incorrect
        
        //dispatch simplified meshlets when simplification is precise enough
        if (groupSimplificationIsPreciseEnough(world_space_bounding_sphere))
        {
            for (uint simplified_meshlet_group_index = 0; simplified_meshlet_group_index < GROUP_SPLIT_COUNT; simplified_meshlet_group_index++)
            {
                queueMeshletForDispatch(current_group.simplified_meshlets[simplified_meshlet_group_index], current_task.scene_object_id, 0); // LOD morphing not yet implemented !!!!!!!!!
            }
        }
        // dispatch base meshlets when current group is a leaf node of the hierarchy tree
        else if (!(current_group.childCount))
        {
            for (uint base_meshlet_group_index = 0; base_meshlet_group_index < current_group.meshlet_count; base_meshlet_group_index++)
            {
                queueMeshletForDispatch(current_group.meshlets[base_meshlet_group_index], current_task.scene_object_id, 0); // LOD morphing not yet implemented !!!!!!!!!
            }
        }
        // append child group into work queue if current simplification isnt precise enough and current group isnt a leaf node
        else
        {
            for (uint child_group_index = 0; child_group_index < current_group.childCount; child_group_index++)
            {
                S_WorkQueueEntry new_task;
                new_task.scene_object_id = current_task.scene_object_id;
                new_task.group_id = current_group.children[child_group_index];
                appendTask(new_task);  
            }
        }
    }
    
    //after work queue is empty dispatch all collected meshlets
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
}
