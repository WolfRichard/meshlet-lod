#pragma once
#include <DirectXMath.h>

using uint		 = unsigned int;
using uint2      = DirectX::XMUINT2;
using uint3      = DirectX::XMUINT3;
using uint4      = DirectX::XMUINT4;
                 
using float2	 = DirectX::XMFLOAT2;
using float3	 = DirectX::XMFLOAT3;
using float4	 = DirectX::XMFLOAT4;
using float4x4	 = DirectX::XMMATRIX;

using quaternion = DirectX::XMVECTOR;

template<typename T>
T lerp(T a, T b, T t) {
    return a + t * (b - a);
}
