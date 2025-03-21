#include "structures.fxh"


// Objects Buffer
StructuredBuffer<SceneObject> objectsBuffer         : register(t0, space0);
// Per Mesh: Meshlet Counts Buffer
StructuredBuffer<uint> meshletCountsBuffer          : register(t1, space0);
// Mesh LoD Structure
StructuredBuffer<MeshLoDStructure> meshLoDStructure : register(t3, space0);
// Draw task arguments
StructuredBuffer<DrawTask> drawTaskBuffers[]        : register(t0, space4);

ConstantBuffer<Constants> constantsBuffer           : register(b0, space0);


cbuffer InstanceIDBuffer                            : register(b1, space0)
{
    uint object_id;
};

cbuffer ObjectLoDBuffer                             : register(b2, space0)
{
    float object_LoD;
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
    float3 directionConeCamera = normalize(constantsBuffer.CameraWorldPos - coneApex);
    return (dot(directionConeCamera, coneAxis) >= coneCutoff);
}

float CalcDetailLevel(float3 center, float radius)
{
    float3 pos = mul(float4(center, 1.0), constantsBuffer.ViewMat).xyz;
    float dist2 = dot(pos, pos);
    
    // calculate the sphere size in screen space
    float size = 2 * constantsBuffer.CoTanHalfFoV * radius / sqrt(dist2 - radius * radius);
    

    float level = clamp(1.0 - size, 0.0, 0.999);
    return level * level;
}


[numthreads(GROUP_SIZE, 1, 1)]
void main(in uint I : SV_GroupIndex,
          in uint wg : SV_GroupID)
{
    // Reset the counter from the first thread in the group
    if (I == 0)
    {
        s_TaskCount = 0;
        s_Payload.object_id = object_id;
    }

    GroupMemoryBarrierWithGroupSync();

    const uint gid = wg * GROUP_SIZE + I;
    
    
    
    //uint object_ID = 0; // incrementingConstant;
    
    uint mesh_lod_index = meshLoDStructure[objectsBuffer[object_id].mesh_id].mesh_offset + object_LoD * meshLoDStructure[objectsBuffer[object_id].mesh_id].lod_count;
    
    
    if (gid < meshletCountsBuffer[mesh_lod_index])
    {
        
        DrawTask task = drawTaskBuffers[mesh_lod_index][gid];
        
        float3 bounding_sphere_center = mul(float4(task.culling_info.bounding_sphere_center, 1.0), objectsBuffer[object_id].object_matrix).xyz;
        float bounding_sphere_radius = length(mul(float4(task.culling_info.bounding_sphere_center + float3(task.culling_info.bounding_sphere_radius, 0, 0), 1.0), objectsBuffer[object_id].object_matrix).xyz - bounding_sphere_center);
        
        float3 normal_cone_apex = mul(float4(task.culling_info.normal_cone_apex, 1.0), objectsBuffer[object_id].object_matrix).xyz;
        float3 normal_cone_axis = mul(float4(task.culling_info.normal_cone_axis, 0.0), objectsBuffer[object_id].object_matrix).xyz;
        
        
        if ((!(constantsBuffer.BoolConstants & FRUSTUM_CULLING_BIT_POS) || IsInFrustum(bounding_sphere_center, bounding_sphere_radius)) &&
            (!(constantsBuffer.BoolConstants & CONE_CULLING_BIT_POS) || !IsShowingBackside(normal_cone_apex, normal_cone_axis, task.culling_info.normal_cone_cutoff)))
        {
        // Acquire an index that will be used to safely access the payload.
        // Each thread gets a unique index.
            uint index = 0;
            InterlockedAdd(s_TaskCount, 1, index);
            
            
            uint meshlet_based_lod = 0;
            if (constantsBuffer.BoolConstants & ENABLE_MESHLET_LOD)
            {
                meshlet_based_lod = CalcDetailLevel(bounding_sphere_center, bounding_sphere_radius) * MESHLET_LOD_COUNT;
            }
            
    
            s_Payload.vertex_count[index] = task.vertex_count[meshlet_based_lod];
            s_Payload.triangle_count[index] = task.triangle_count[meshlet_based_lod];
    
            s_Payload.vertex_offset[index] = task.vertex_offset[meshlet_based_lod];
            s_Payload.triangle_offset[index] = task.triangle_offset[meshlet_based_lod];
            
            s_Payload.meshlet_LoD[index] = meshlet_based_lod;
            
            
        }
    
        
    }
   
    // All threads must complete their work so that we can read s_TaskCount
    GroupMemoryBarrierWithGroupSync();
    DispatchMesh(s_TaskCount, 1, 1, s_Payload);
}
