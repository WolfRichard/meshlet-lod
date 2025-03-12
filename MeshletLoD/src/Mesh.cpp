#include "Mesh.h"
#include "Animator.h"
#include <deque>

using namespace DirectX;

Mesh::Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
    
    
    parseMesh(assimp_mesh, assimp_scene);

}


void Mesh::parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene)
{
   

    
}

void Mesh::generateMeshlets()
{
    
}
