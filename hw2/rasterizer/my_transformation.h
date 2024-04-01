//
// Created by Bengisu on 12/12/2023.
//

#ifndef CODE_TEMPLATE_MY_TRANSFORMATION_H
#define CODE_TEMPLATE_MY_TRANSFORMATION_H

#include "Scene.h"
#include "Mesh.h"
#include "Matrix4.h"
#include "Helpers.h"
#include "my_line.h"


Matrix4 modeling_transformations(const Scene& scene, Mesh& mesh);
void culling(const Camera& cam, Mesh& mesh, vector<int>& visible_triangles, bool cullingEnabled);


#endif //CODE_TEMPLATE_MY_TRANSFORMATION_H
