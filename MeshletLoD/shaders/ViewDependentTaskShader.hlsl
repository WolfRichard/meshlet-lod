#include "ViewDependentStructures.fxh"


ConstantBuffer<S_Constants> constants                       : register(b0, space0);
StructuredBuffer<S_SceneObject> objectsBuffer               : register(t0, space0);

RWStructuredBuffer<S_WorkQueueEntry> workQueue              : register(u0, space0);
RWStructuredBuffer<S_WorkQueueCounters> workQueueCounters   : register(u1, space0);
RWStructuredBuffer<S_PayloadEntry> payloadBuffer            : register(u2, space0);


// bindless buffers
StructuredBuffer<S_Meshlet> meshletBuffers[]                : register(t0, space1);


// Payload will be used in the mesh shader.
groupshared S_Payload gs_Payload;

// The number of meshlets that are dispatched by this thread group,
groupshared uint gs_MeshletCount;


void appendTask(S_WorkQueueEntry new_task)
{
    uint task_index = 0;
    InterlockedAdd(workQueueCounters[0].work_queue_tail, 1, task_index);    // Atomically increment endCounter
    workQueue[task_index % WORK_QUEUE_SIZE] = new_task;     // loop index because of ring buffer structure
}


bool consumeTask(out S_WorkQueueEntry out_task)
{
    uint task_index, prev_index;
    do
    {
        prev_index = workQueueCounters[0].work_queue_head;
        uint end = workQueueCounters[0].work_queue_tail;

        if (prev_index >= end)
            return false; // No work available
        
        // Try to atomically update 'beginCounter' only if it hasn’t changed
        InterlockedCompareExchange(workQueueCounters[0].work_queue_head, prev_index, prev_index + 1, task_index);
        
    } while ((task_index != prev_index)); // lost race condition in between checking for available work and incrementing counter --> so check if there is other work available to grab
    
    // Load the task safely as thread won race condition after checking that there was still work to do
    out_task = workQueue[task_index % WORK_QUEUE_SIZE];
    return true;
}

// Frustum culling in world space
bool isInFrustum(float3 center, float radius)
{
    if (!(constants.BoolConstants & FRUSTUM_CULLING_BIT_POS))
        return true;
    // TODO: FIX CULLING currently cuts away to early
    float4 f4Center = float4(center, 1.0);
    for (int i = 0; i < 6; ++i)
    {
        if (dot(constants.Frustum[i], f4Center) < -radius - 1)
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

bool groupSimplificationIsPreciseEnough(S_BoundingSphere bounding_sphere, uint lod_level, out uint tessellation_level) // bounding sphere must be in world space!
{
    float cam_dist = max(distance(constants.CameraWorldPos, bounding_sphere.center) - bounding_sphere.radius, 0.000000001);
    float expected_lod = log2(cam_dist / constants.LoD_Scale);
    tessellation_level = uint(min(max(expected_lod * -1 + 1, 0), MAX_TESSELLATION_LEVEL));
    return lod_level <= max(expected_lod, 0);
}


bool isPreciseEnough(S_BoundingSphere bounding_sphere) // bounding sphere must be in clip space
{
    float d2 = dot(bounding_sphere.center, bounding_sphere.center);
    float r2 = bounding_sphere.radius * bounding_sphere.radius;
    float sphere_diameter_uv = max(constants.ProjMat[0][0], constants.ProjMat[1][1]) * bounding_sphere.radius / sqrt(d2 - r2);
    float view_size = max(constants.ScreenWidth, constants.ScreenHeight);
    return (sphere_diameter_uv * view_size) < constants.LoD_Scale; //    1.0;
}


void queueMeshletForDispatch(uint meshlet_index, uint object_index, uint tessellation_level, uint meshlet_triangle_count)
{
    if (tessellation_level == 0)
    {
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

void processTask(S_WorkQueueEntry task)
{
    S_SceneObject scene_object = objectsBuffer[task.scene_object_id];
    S_Meshlet current_meshlet = meshletBuffers[scene_object.mesh_id][task.meshlet_id];
        
    bool precise_enough;
            
    S_BoundingSphere world_space_bounding_sphere;
    world_space_bounding_sphere.center = mul(float4(current_meshlet.bounding_sphere.center, 1.0), scene_object.object_matrix).xyz; // object --> world space
    
    uint tessellation_level = 0;
        
    if (constants.BoolConstants & SCREEN_SPACE_ERROR_BASED_LOD_BIT_POS)
    {
        world_space_bounding_sphere.radius = current_meshlet.base_error * ExtractMaxScaleFactor(scene_object.object_matrix);
        world_space_bounding_sphere.center = mul(float4(world_space_bounding_sphere.center, 1), constants.ViewMat).xyz; // world --> clip space
        precise_enough = isPreciseEnough(world_space_bounding_sphere); // actually in screen clip space
    }
    else
        world_space_bounding_sphere.radius = current_meshlet.bounding_sphere.radius * ExtractMaxScaleFactor(scene_object.object_matrix);
        precise_enough = groupSimplificationIsPreciseEnough(world_space_bounding_sphere, current_meshlet.discrete_level_of_detail, tessellation_level);
        
    if (precise_enough)
    {
        if (isInFrustum(world_space_bounding_sphere.center, world_space_bounding_sphere.radius))
            queueMeshletForDispatch(task.meshlet_id, task.scene_object_id, tessellation_level, current_meshlet.triangle_count);
    }
    else // append all child meshlets into work queue if current simplification isnt precise enough (zero children when current meshlet is a leaf node)
        for (uint child_meshlet_index_index = 0; child_meshlet_index_index < current_meshlet.child_count; child_meshlet_index_index++)
        {
            S_WorkQueueEntry new_task;
            new_task.scene_object_id = task.scene_object_id;
            new_task.meshlet_id = current_meshlet.child_meshlets[child_meshlet_index_index];
            appendTask(new_task);
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
       
    // processing of work queue entries    
    S_WorkQueueEntry current_task;
    while (consumeTask(current_task))
    {
        processTask(current_task);
    }
      
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
