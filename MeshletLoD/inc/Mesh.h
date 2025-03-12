#pragma once

#include <meshoptimizer.h>

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

#include "HLSLnames.h"

#include <chrono>
#include <map>

namespace
{

#include "../shaders/ViewDependentStructures.fxh"

} // namespace

class Mesh
{
public:
    Mesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);



private:
    void parseMesh(aiMesh* assimp_mesh, const aiScene* assimp_scene);
    void generateMeshlets();
    
};
