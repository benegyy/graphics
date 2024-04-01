#include <iostream>
#include <vector>
#include <string>
#include "Scene.h"
#include "Matrix4.h"
#include "Helpers.h"
#include "my_transformation.h"
#include "Triangle.h"
#include "Scene.h"

using namespace std;

Scene *scene;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Please run the rasterizer as:" << endl
             << "\t./rasterizer <input_file_name>" << endl;
        return 1;
    }
    else
    {
        const char *xmlPath = argv[1];

        scene = new Scene(xmlPath);
        for(Mesh* mesh: scene->meshes)
        {
            Matrix4 mModel = modeling_transformations(*scene, *mesh);
            mesh->transformedVertices = {};
            for(Triangle triangle: mesh->triangles)
            {
                Vec4 v1,v2,v3;

                v1 = convertVec4(*(scene->vertices[triangle.vertexIds[0]-1]));
                v2 = convertVec4(*(scene->vertices[triangle.vertexIds[1]-1]));
                v3 = convertVec4(*(scene->vertices[triangle.vertexIds[2]-1]));

                mesh->transformedVertices[triangle.vertexIds[0]-1] = multiplyMatrixWithVec4(mModel, v1);
                mesh->transformedVertices[triangle.vertexIds[1]-1] = multiplyMatrixWithVec4(mModel, v2);
                mesh->transformedVertices[triangle.vertexIds[2]-1] = multiplyMatrixWithVec4(mModel, v3);


            }
        }
        for (int i = 0; i < scene->cameras.size(); i++)
        {
            // initialize image with basic values
            scene->initializeImage(scene->cameras[i]);

            // do forward rendering pipeline operations
            scene->forwardRenderingPipeline(scene->cameras[i]);

            // generate PPM file
            scene->writeImageToPPMFile(scene->cameras[i]);

            // Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
            // Change/remove implementation if necessary.
            scene->convertPPMToPNG(scene->cameras[i]->outputFilename);
        }

        return 0;
    }
}