// Input buffers (read-only)
StructuredBuffer<float> InputA : register(t0);
StructuredBuffer<float> InputB : register(t1);

// Output buffer (read-write)
RWStructuredBuffer<float> Output : register(u0);

// Compute shader entry point
[numthreads(256, 1, 1)] // Defines the number of threads per thread group
void main(uint3 dispatchThreadID : SV_DispatchThreadID)
{
    // Fetch the thread ID
    uint index = dispatchThreadID.x;

    // Perform element-wise addition
    Output[index] = InputA[index] + InputB[index];
}