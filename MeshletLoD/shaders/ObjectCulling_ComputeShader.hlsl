# include "Structures.fxh"

// Input buffers (read-only)

//Constants Buffer for frustum information
cbuffer InstanceIDBuffer : register(b0, space0)
{
    Constants constantsBuffer;
};

//Constants Buffer for frustum information
cbuffer InstanceIDBuffer : register(b1, space0)
{
    uint maxObjectCount;
};

// Objects Buffer
StructuredBuffer<SceneObject> objectsBuffer : register(t0, space0);
// Per Mesh: Meshlet Counts Buffer
StructuredBuffer<uint> meshletCountsBuffer : register(t1, space0);


// Output buffer (read-write)

// number of visible objects
RWStructuredBuffer<uint> visibleObjectCount : register(u0, space0);

// arguments buffer for dispatch indirect
RWStructuredBuffer<CommandStructure> indirectArgumentBuffer : register(u1, space0);



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


float CalcDetailLevel(float3 center, float radius)
{
    // vector from camera to object center
    float3 pos = mul(float4(center, 1.0), constantsBuffer.ViewMat).xyz;
    
    // Square of distance from camera to circumscribed sphere
    float dist2 = dot(pos, pos);
    
    // Calculate the sphere size in screen space
    float size = constantsBuffer.CoTanHalfFoV * radius / sqrt(dist2 - radius * radius);
    
    // Calculate detail level
    float level = clamp(1.0 - size, 0.0, 1.0);
    return level;
}


[numthreads(GROUP_SIZE, 1, 1)]
void main(uint3 dispatchThreadID : SV_DispatchThreadID)
{    
    uint objectID = dispatchThreadID.x;
 
    // idle all threads without workload
    if (objectID >= maxObjectCount)
        return;
    
    // Reset the counter from the first thread in the group
    if (objectID == 0)
    {
        visibleObjectCount[0] = 0;
    }
    GroupMemoryBarrierWithGroupSync();
    
    
    SceneObject object = objectsBuffer[objectID];
    if (IsInFrustum(object.bounding_sphere_center, object.bounding_sphere_radius))
    {
        uint argumentBufferID;
        InterlockedAdd(visibleObjectCount[0], 1, argumentBufferID);
        uint thread_count = ((meshletCountsBuffer[object.mesh_id] + GROUP_SIZE - 1) / GROUP_SIZE) * GROUP_SIZE;
        
        indirectArgumentBuffer[argumentBufferID].instanceID = objectID;
        indirectArgumentBuffer[argumentBufferID].level_of_detail = CalcDetailLevel(object.bounding_sphere_center, object.bounding_sphere_radius);
        indirectArgumentBuffer[argumentBufferID].dispatchArguments = uint3(thread_count, 1, 1);       
    }
}