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
    uint task_index = 0;
    InterlockedAdd(workQueueCounters[1], 1, task_index);    // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task;     // loop index because of ring buffer structure
}


// TODO: POSIBLE RACE CONDITION! should swap to InterlockedCompareExchange() instead!!!!!!!!!!!!!!!!!!!!!!!!!1
bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint task_index, prev_index;
    do
    {
        prev_index = workQueueCounters[0]; // Read beginCounter
        uint end = workQueueCounters[1]; // Read endCounter

        if (prev_index >= end)
            return false; // No work available
        
        // Try to atomically update 'beginCounter' only if it hasn’t changed
        InterlockedCompareExchange(workQueueCounters[0], prev_index, prev_index + 1, task_index);
        
    } while (task_index != prev_index); // lost race condition in between checking for available work and incrementing counter --> so check if there is other work available to grab
    
    
    // Load the task safely as thread won race condition after checking that there was still work to do
    out_task = workQueue[task_index % WORK_QUEUE_SIZE];
    return true;
}


// Payload will be used in the mesh shader.
groupshared S_Payload gs_Payload;

// The number of meshlets that are dispatched by this thread group,
groupshared uint gs_MeshletCount;


// Frustum culling in world space
bool isInFrustum(float3 center, float radius)
{
    return true;
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


bool groupSimplificationIsPreciseEnough(S_BoundingSphere bounding_sphere, float lod_error) // bounding sphere must be in world space!
{
    float cam_dist = max(distance(constants.CameraWorldPos, bounding_sphere.center) - bounding_sphere.radius, 0);
    //float cam_dist = distance(constants.CameraWorldPos, bounding_sphere.center);
    float lod_threshold = max(log2(cam_dist / constants.LoD_Scale), 0);
    return lod_error <= lod_threshold;
}

bool errorLessThanPixel(S_BoundingSphere bounding_sphere) // bounding sphere must be in clip space
{
    float d2 = dot(bounding_sphere.center, bounding_sphere.center);
    float r2 = bounding_sphere.radius * bounding_sphere.radius;
    float sphere_diameter_uv = max(constants.ProjMat[0][0], constants.ProjMat[1][1]) * bounding_sphere.radius / sqrt(d2 - r2);
    float view_size = max(constants.ScreenWidth, constants.ScreenHeight);
    float sphere_diameter_pixels = sphere_diameter_uv * view_size;
    return sphere_diameter_pixels < 1.0;
}

void queueMeshletForDispatch(uint meshlet_index, uint object_index, float lod_blend_value, uint lod_depth)
{
    S_SceneObject scene_object = objectsBuffer[object_index];
    S_BoundingSphere world_space_bounding_sphere;
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][meshlet_index];
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz;
    world_space_bounding_sphere.radius = current_meshlet.bounding_sphere.radius * ExtractMaxScaleFactor(scene_object.object_matrix);
    if (isInFrustum(world_space_bounding_sphere.center, world_space_bounding_sphere.radius))
    {
        uint payload_index = 0;
        InterlockedAdd(gs_MeshletCount, 1, payload_index); // Undefined behavior if max meshlet count per work group is set to low 
        gs_Payload.tasks[payload_index].meshlet_id = meshlet_index;
        gs_Payload.tasks[payload_index].object_id = object_index;
        gs_Payload.tasks[payload_index].lod_morphing = lod_blend_value;
        gs_Payload.tasks[payload_index].lod_tree_depth = lod_depth;
    }
}

[numthreads(GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
    uint global_thread_index = I + wg * GROUP_SIZE;
   
    
    // debug (render most simplified meshlets)
    /*
    if (I == 0)
    {
        gs_MeshletCount = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    if (global_thread_index == 0)
    {
        S_SceneObject s_o = objectsBuffer[0];
        S_MeshletGroup root_group = groupBuffers[s_o.mesh_id][s_o.root_group_id];
        queueMeshletForDispatch(root_group.simplified_meshlets[0], 0, 0);
        queueMeshletForDispatch(root_group.simplified_meshlets[1], 0, 0);
    }
  
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
    return;
    */

    
    
   // Reset the group shared counter variables from the first thread in the group
    if (I == 0)
    {
        gs_MeshletCount = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    
    /*
    // object culling and extracting of root nodes into work queue
    for (uint scene_object_index = global_thread_index; scene_object_index < constants.SceneObjectCount; scene_object_index += PERSISTENT_THREAD_COUNT)
    {
        S_SceneObject scene_object = objectsBuffer[scene_object_index];
        
        if (isInFrustum(scene_object.bounding_sphere.center, scene_object.bounding_sphere.radius))
        {
            S_WorkQueueEntry new_task;
            new_task.group_id = scene_object.root_group_id;
            new_task.scene_object_id = scene_object_index;
            appendTask(new_task);
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
        if (groupSimplificationIsPreciseEnough(world_space_bounding_sphere, current_group.hierarchy_tree_depth))
        {
            for (uint simplified_meshlet_group_index = 0; simplified_meshlet_group_index < GROUP_SPLIT_COUNT; simplified_meshlet_group_index++)
            {
                queueMeshletForDispatch(current_group.simplified_meshlets[simplified_meshlet_group_index], current_task.scene_object_id, current_task.group_id, current_group.hierarchy_tree_depth + 1); // TODO: LOD morphing not yet implemented currently used for debug !!!!!!!!!
            }
        }
        // dispatch base meshlets when current group is a leaf node of the hierarchy tree
        else if (!(current_group.childCount))
        {
            for (uint base_meshlet_group_index = 0; base_meshlet_group_index < current_group.meshlet_count; base_meshlet_group_index++)
            {
                queueMeshletForDispatch(current_group.meshlets[base_meshlet_group_index], current_task.scene_object_id, current_task.group_id, current_group.hierarchy_tree_depth); // TODO:  LOD morphing not yet implemented currently used for debug!!!!!!!!!
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
    */
    
    // debug test (rendering first object without tree traversal by evaluating every meshlet simultaniously)
    S_SceneObject scene_object = objectsBuffer[0];
    if (global_thread_index < scene_object.mesh_meshlet_count)
    {
        float world_scale = ExtractMaxScaleFactor(scene_object.object_matrix);
        
        S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][global_thread_index];
        S_BoundingSphere base_bounding_sphere;
        base_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz; // object --> world space
        base_bounding_sphere.radius = current_meshlet.base_error * world_scale;
        base_bounding_sphere.center = mul(float4(base_bounding_sphere.center, 1), constants.ViewMat).xyz; // world --> clip space
        
        
        S_BoundingSphere simplified_bounding_sphere;
        simplified_bounding_sphere.center = mul(float4(current_meshlet.simplified_group_bounds.center, 1.0), scene_object.object_matrix).xyz; // object --> world space
        simplified_bounding_sphere.radius = current_meshlet.simplification_error * world_scale;
        simplified_bounding_sphere.center = mul(float4(simplified_bounding_sphere.center, 1), constants.ViewMat).xyz; // world --> clip space
        
        bool parent_precise_enough = errorLessThanPixel(simplified_bounding_sphere);
        bool base_precise_enough = errorLessThanPixel(base_bounding_sphere);
        

        if (!parent_precise_enough && base_precise_enough)
        {
            queueMeshletForDispatch(global_thread_index, 0, current_meshlet.group_id, current_meshlet.discrete_level_of_detail);
        }
    }
    
    
    
    
    
    
    //after work queue is empty dispatch all collected meshlets
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
}
