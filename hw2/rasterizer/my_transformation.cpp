#include "my_transformation.h"
#include "Helpers.h"
#include "my_line.h"

#include <vector>

using namespace std;

Matrix4 modeling_transformations(const Scene& scene, Mesh& mesh)
{
    Matrix4 result = getIdentityMatrix();
    Matrix4 temp;
    vector<int>::iterator id = mesh.transformationIds.begin();
    vector<char>::iterator type = mesh.transformationTypes.begin();
    for( ; id != mesh.transformationIds.end(); id++, type++){
        int transfromationID=*id-1;
        switch(*type)
        {
            case 'r':
                temp = scene.rotations[transfromationID]->rotate_matrix();
                break;
            case 't':
                temp = scene.translations[transfromationID]->translate_matrix();
                break;
            case 's':
                temp = scene.scalings[transfromationID]->scale_matrix();
                break;
        }
        result = multiplyMatrixWithMatrix(temp, result);

    }
    return result;
    }


void culling(const Camera& cam, Mesh& mesh, vector<int>& visible_triangles, bool cullingEnabled)
{
    if(cullingEnabled)
    {
        size_t i=0;
        for(Triangle triangle: mesh.triangles)
        {
            Vec3 corner1,corner2,corner3;
            Vec3 edge1, edge2;
            Vec3 triangle_normal;// triangle_normal of the triangle
            Vec3 pos = cam.position;

            corner1 = convertVec3(mesh.transformedVertices[triangle.vertexIds[0] - 1]);
            corner2 = convertVec3(mesh.transformedVertices[triangle.vertexIds[1] - 1]);
            corner3 = convertVec3(mesh.transformedVertices[triangle.vertexIds[2] - 1]);

            edge1 = subtractVec3(corner2, corner1);
            edge2 = subtractVec3(corner3, corner2);

            triangle_normal = crossProductVec3(edge1, edge2);

            Vec3 vStare = subtractVec3(corner1 , pos);

            double dotProduct = dotProductVec3(triangle_normal, vStare);
            if(dotProduct  < 0)
            {
                visible_triangles.push_back(i);
            }

            i++;
        }
    }
    else
    {
        for(size_t i=0; i<mesh.triangles.size(); i++)
        {
            visible_triangles.push_back(i);
        }
    }
}
