#include "Structures.fxh"


// read only buffers
ConstantBuffer<S_Constants> constants                       : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer               : register(t0, space0);

// read-write buffers
RWStructuredBuffer<S_WorkQueueEntry> workQueue              : register(u0, space0);
RWStructuredBuffer<S_WorkQueueCounters> workQueueCounters   : register(u1, space0);
RWStructuredBuffer<S_PayloadEntry> payloadBuffer            : register(u2, space0);

// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]                : register(t0, space1);

// group shared variables
groupshared S_Payload gs_Payload; // payload that will be passed to the mesh shader.
groupshared uint gs_MeshletCount; // The number of meshlets that are dispatched by this thread group.



// pushes the given task into the FIFO work queue
// and atomically adjusts the tail pointer accordingly
void appendTask(S_WorkQueueEntry new_task)
{
    uint task_index = 0;
    InterlockedAdd(workQueueCounters[0].work_queue_tail, 1, task_index);    // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task;                     // loop index because of ring buffer structure
}


// tries to retrieve a task from the work queue, while returning false if the queue was empty
// to save resources only the first lane of a warp allocates upt to GROUP_SIZE tasks which are then distributed to the other threads
// it also updates the head pointer of the work queue accordingly (also atomic)
bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint group_assigned_work_count = 0;     // counter of how much tasks the first lane of the wave did manage to claim for the thread group
    uint work_queue_task_group_offset = 0;  // index to the first task inside the work queue that got claimed by the group
    
    if (WaveIsFirstLane()) // only first thread of a group tries to catch work, and then dsitributes it to the whole group
    {
        uint task_index, prev_index;
        do
        {
            prev_index = workQueueCounters[0].work_queue_head;
            uint end = workQueueCounters[0].work_queue_tail;

            if (prev_index >= end)
                break; // No work available, break out of the while loop
        
            // only claim up to one task per thread of the work group
            uint intended_work_claim_count = min(end - prev_index, GROUP_SIZE);
        
            // Try to atomically update the head pointer, but only if it hasn’t changed
            InterlockedCompareExchange(workQueueCounters[0].work_queue_head, prev_index, prev_index + intended_work_claim_count, task_index);
        
            group_assigned_work_count = (task_index == prev_index) ? intended_work_claim_count : 0;
        
        } while (group_assigned_work_count == 0); // lost race condition in between checking for available work and incrementing counter --> so check if there is other work available to grab
        
        work_queue_task_group_offset = task_index;
    }
    
    // broadcast tasks from first lane to the whole wave
    group_assigned_work_count = WaveReadLaneFirst(group_assigned_work_count);
    work_queue_task_group_offset = WaveReadLaneFirst(work_queue_task_group_offset);

    if (WaveGetLaneIndex() >= group_assigned_work_count)
        return false;
    
    // Load the task safely as thread won race condition after checking that there was still work to do
    out_task = workQueue[(work_queue_task_group_offset + WaveGetLaneIndex()) % WORK_QUEUE_SIZE];
    return true;
}


// perform frustum culling check on a sphere in world space
// frustum is defined as 6 planes within the constant buffer
bool isInFrustum(float3 center, float radius)
{
    if (!(constants.BoolConstants & FRUSTUM_CULLING_BIT_POS))
        return true;
    // TODO: FIX CULLING currently cuts away to early
    float4 f4Center = float4(center, 1.0);
    for (int i = 0; i < 6; ++i)
    {
        if (dot(constants.Frustum[i], f4Center) < -radius - 1) // TODO: the "-1" is just a bandage-fix
            return false;
    }
    return true;
}


// returns the highest scale factor considering each axis of a matrix
// usefull for conservative radius estimation when transforming a sphere into another coordinate space 
float ExtractMaxScaleFactor(float4x4 m)
{
    float sx = length(m[0].xyz);
    float sy = length(m[1].xyz);
    float sz = length(m[2].xyz);
   
    return max(sx, max(sy, sz));
}


// returns if the provided bounding sphere would be precise enough for the given camera position
// bounding sphere parameter has to be passed in world space coordiantes
// also returns the tessellation level that would be necessary for the given sphere
bool groupSimplificationIsPreciseEnough(S_BoundingSphere bounding_sphere, uint lod_level, out uint tessellation_level)
{
    //return lod_level < constants.DebugFloatSliderValue;
    
    float cam_dist = max(distance(constants.CameraWorldPos, bounding_sphere.center) - bounding_sphere.radius, 0.000000001);
    float expected_lod = log2(cam_dist / constants.LoD_Scale);
    tessellation_level = uint(min(max(expected_lod * -1 + 1, 0), MAX_TESSELLATION_LEVEL));
    return lod_level <= max(expected_lod, 0);
}



// similar to groupSimplificationIsPreciseEnough() but uses screen space size instead of camera distance
// the provided bounding sphere has to be provided in clip space
bool isPreciseEnough(S_BoundingSphere bounding_sphere)
{
    float d2 = dot(bounding_sphere.center, bounding_sphere.center);
    float r2 = bounding_sphere.radius * bounding_sphere.radius;
    float sphere_diameter_uv = max(constants.ProjMat[0][0], constants.ProjMat[1][1]) * bounding_sphere.radius / sqrt(d2 - r2);
    float view_size = max(constants.ScreenWidth, constants.ScreenHeight);
    // through the LoD_Scale constant, screen space error is allowed to be bigger than 1 pxl but this can result in visual pop-in artefacts
    return (sphere_diameter_uv * view_size) < constants.LoD_Scale; 
}


// pushes a certain number of payload entries into the global payload buffer
// the number is automatically adjusted according to tessellation grade, such that the limitations of the mesh shader output is not exceeded when generating aditional primitives
void queueMeshletForDispatch(uint meshlet_index, uint object_index, uint tessellation_level, uint meshlet_triangle_count)
{
    if (tessellation_level == 0)
    {
        // as there is no tessellation planed for this meshlet, the original structure can be directly used
        uint local_meshlet_count = 0;
        InterlockedAdd(gs_MeshletCount, 1, local_meshlet_count);
    
        uint global_payload_index = 0;
        InterlockedAdd(workQueueCounters[0].payload_tail, 1, global_payload_index);
        payloadBuffer[global_payload_index].meshlet_id = meshlet_index;
        payloadBuffer[global_payload_index].object_id = object_index;
        payloadBuffer[global_payload_index].tessellation_grade = 0;
        payloadBuffer[global_payload_index].tessellation_triangle_offset = 0;
        payloadBuffer[global_payload_index].tessellation_triangle_count = 0;
    }
    else
    {
        // for tessellated meshlets, use the look up table to check into how many sub-meshlet-tasks it has to be divided
        uint triangles_per_instance = TESSELLATION_INITIAL_TRIANGLE_COUNTS[tessellation_level - 1];
        uint required_mesh_shader_instances = (meshlet_triangle_count + triangles_per_instance - 1) / triangles_per_instance;
        
        uint local_meshlet_count = 0;
        InterlockedAdd(gs_MeshletCount, required_mesh_shader_instances, local_meshlet_count);

        uint global_payload_index = 0;
        InterlockedAdd(workQueueCounters[0].payload_tail, required_mesh_shader_instances, global_payload_index);
       
        
        for (uint i = 0; i < required_mesh_shader_instances - 1; i++)
        {
            payloadBuffer[global_payload_index + i].meshlet_id = meshlet_index;
            payloadBuffer[global_payload_index + i].object_id = object_index;
            payloadBuffer[global_payload_index + i].tessellation_grade = tessellation_level;
            payloadBuffer[global_payload_index + i].tessellation_triangle_offset = triangles_per_instance * i;
            payloadBuffer[global_payload_index + i].tessellation_triangle_count = triangles_per_instance;
            
        }
        payloadBuffer[global_payload_index + required_mesh_shader_instances - 1].meshlet_id = meshlet_index;
        payloadBuffer[global_payload_index + required_mesh_shader_instances - 1].object_id = object_index;
        payloadBuffer[global_payload_index + required_mesh_shader_instances - 1].tessellation_grade = tessellation_level;
        payloadBuffer[global_payload_index + required_mesh_shader_instances - 1].tessellation_triangle_offset = triangles_per_instance * (required_mesh_shader_instances - 1);
        payloadBuffer[global_payload_index + required_mesh_shader_instances - 1].tessellation_triangle_count = meshlet_triangle_count % triangles_per_instance;
    }
}


// evaluate if the currently worked on task is precise enough to be rendered to the screen, if not add its children back into the work queue
void processTask(S_WorkQueueEntry task)
{
    S_SceneObject scene_object = objectsBuffer[task.scene_object_id];
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][task.meshlet_id];
        
    bool precise_enough;
            
    S_BoundingSphere world_space_bounding_sphere;
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz; // object --> world space
    world_space_bounding_sphere.radius = current_meshlet.bounding_sphere.radius * ExtractMaxScaleFactor(scene_object.object_matrix);
    
    uint tessellation_level = 0;
        
    if (constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS)
    {
        S_BoundingSphere clip_space_bounding_sphere;
        clip_space_bounding_sphere.radius = current_meshlet.base_error * ExtractMaxScaleFactor(scene_object.object_matrix);
        clip_space_bounding_sphere.center = mul(float4(world_space_bounding_sphere.center, 1), constants.ViewMat).xyz; // world --> clip space
        precise_enough = isPreciseEnough(clip_space_bounding_sphere); 
    }
    else
        precise_enough = groupSimplificationIsPreciseEnough(world_space_bounding_sphere, current_meshlet.discrete_level_of_detail, tessellation_level);
        
    if (precise_enough)
    {
        if (isInFrustum(world_space_bounding_sphere.center, world_space_bounding_sphere.radius))
            queueMeshletForDispatch(task.meshlet_id, task.scene_object_id, tessellation_level, current_meshlet.triangle_count);
    }
    else // append all child meshlets into work queue if current simplification isnt precise enough
        for (uint child_meshlet_index_index = 0; child_meshlet_index_index < current_meshlet.child_count; child_meshlet_index_index++)
        {
            S_WorkQueueEntry new_task;
            new_task.scene_object_id = task.scene_object_id;
            new_task.meshlet_id = current_meshlet.child_meshlets[child_meshlet_index_index];
            appendTask(new_task);
        }
}


// task shader entry point
[shader("amplification")]
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
       
    // processing of work queue entries    
    S_WorkQueueEntry current_task;
    bool got_work;
    do
    {
        got_work = consumeTask(current_task);
        if (got_work)
        {
            processTask(current_task);
        }
    }
    while (WaveActiveAnyTrue(got_work));
      
    //after work queue is empty dispatch all collected meshlets
    GroupMemoryBarrierWithGroupSync();
    if (I == 0)
    {
        uint global_offset = 0;
        InterlockedAdd(workQueueCounters[0].payload_head, gs_MeshletCount, global_offset);
        gs_Payload.global_payload_offset = global_offset;
    }
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(gs_MeshletCount, 1, 1, gs_Payload);
}
