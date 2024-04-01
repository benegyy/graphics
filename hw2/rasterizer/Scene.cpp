#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>


#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include "my_transformation.h"
#include "my_line.h"
#include <map>
#include <unordered_map>
#include <limits>

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(INFINITY);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = INFINITY;
				this->depth[i][j] = INFINITY;
				this->depth[i][j] = INFINITY;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}
void updateDepthBuffer(const Vec4& vertex, vector<vector<double>>& depthBuffer) {
    double x = static_cast<int>(vertex.x);
    double y = static_cast<int>(vertex.y);
    if (x >= 0 && x < depthBuffer[0].size() && y >= 0 && y < depthBuffer.size()) {
        double depth = vertex.z;
        depthBuffer[y][x] = min(depthBuffer[y][x], depth);
    }
}
bool isPointInsideTriangle(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3) {
    float d1 = (x - x2) * (y1 - y2) - (x1 - x2) * (y - y2);
    float d2 = (x - x3) * (y2 - y3) - (x2 - x3) * (y - y3);
    float d3 = (x - x1) * (y3 - y1) - (x3 - x1) * (y - y1);
    return (d1 >= 0 && d2 >= 0 && d3 >= 0) ;
}

void Scene::forwardRenderingPipeline(Camera *camera) {
 //
    vector<vector<double>> depthBuffer;
    depthBuffer.resize(camera->verRes, vector<double>(camera->horRes, numeric_limits<double>::infinity()));

//
    for (Mesh *mesh: meshes) {
        vector<int> visible_triangle = {};
        unordered_map<int, Vec4> vertices_transformed = {};

        culling(*camera, *mesh, visible_triangle, cullingEnabled);


        Matrix4 MperMcam = camera->computeMperMcam();

        for (int triangleId: visible_triangle) {
            Vec4 v1, v2, v3;
            Triangle triangle = mesh->triangles[triangleId];

            v1 = mesh->transformedVertices.at(triangle.vertexIds[0] - 1);
            v2 = mesh->transformedVertices.at(triangle.vertexIds[1] - 1);
            v3 = mesh->transformedVertices.at(triangle.vertexIds[2] - 1);

            v1 = multiplyMatrixWithVec4(MperMcam, v1);
            v2 = multiplyMatrixWithVec4(MperMcam, v2);
            v3 = multiplyMatrixWithVec4(MperMcam, v3);
//
            if (camera->projectionType) {
                v1 = pers_divide(v1);
                v2 = pers_divide(v2);
                v3 = pers_divide(v3);
            }
            updateDepthBuffer(v1, depthBuffer);
            updateDepthBuffer(v2, depthBuffer);
            updateDepthBuffer(v3, depthBuffer);
            //
            vertices_transformed.emplace(triangle.vertexIds[0], v1);
            vertices_transformed.emplace(triangle.vertexIds[1], v2);
            vertices_transformed.emplace(triangle.vertexIds[2], v3);

        }


        if (!mesh->type)
        {
            vector<my_line> lines = {};
            map<pair<int, int>, bool> isLineAdded = {};

            for (int triangleId: visible_triangle) {
                Triangle triangle = mesh->triangles[triangleId];
                int vertex1, vertex2, vertex3;


                vertex1 = triangle.vertexIds[0];
                vertex2 = triangle.vertexIds[1];
                vertex3 = triangle.vertexIds[2];

                if (vertex2 < vertex1) {
                    swap(vertex1, vertex2);
                }
                if (vertex3 < vertex1) {
                    swap(vertex1, vertex3);
                }
                if (vertex3 < vertex2) {
                    swap(vertex2, vertex3);
                }



                if (!isLineAdded.count(make_pair(vertex1, vertex2))) {
                    Vec4 linep1, linep2;
                    linep1 = vertices_transformed.at(vertex1);
                    linep2 = vertices_transformed.at(vertex2);
                    lines.push_back(my_line(linep1, linep2));
                    isLineAdded.emplace(std::pair<int, int>{vertex1, vertex2}, true);
                }
                if (!isLineAdded.count(make_pair(vertex2, vertex3))) {
                    Vec4 p2, p3;
                    p2 = vertices_transformed.at(vertex2);
                    p3 = vertices_transformed.at(vertex3);

                    lines.push_back(my_line(p2, p3));
                    isLineAdded.emplace(std::make_pair(vertex2, vertex3), true);

                }
                if (!isLineAdded.count(make_pair(vertex1, vertex3))) {
                    Vec4 p1, p3;
                    p1 = vertices_transformed.at(vertex1);
                    p3 = vertices_transformed.at(vertex3);
                    lines.push_back(my_line(p1, p3));
                    isLineAdded.emplace(std::make_pair(vertex1, vertex3), true);
                }
            }

            Matrix4 view_matrix = camera->viewport_transformation();

            for (my_line &line: lines) {
                Vec4 point1, point2;

                if (camera->projectionType) {
                    point1 = pers_divide(line.p1);
                    point2 = pers_divide(line.p2);
                }
                line.p1 = multiplyMatrixWithVec4(view_matrix, point1);
                line.p2 = multiplyMatrixWithVec4(view_matrix, point2);
            }

            for (my_line line: lines) {
                line_algorithm(line, this->image, this->colorsOfVertices);
            }
        } else
        {
            Matrix4 view_matrix = camera->viewport_transformation();

            for (pair<const int, Vec4> &vec: vertices_transformed) {
                if (camera->projectionType) {
                    vec.second = pers_divide(vec.second);
                }
                vec.second = multiplyMatrixWithVec4(view_matrix, vec.second);
            }

            for (int id: visible_triangle) {
                Triangle triangle = mesh->triangles[id];
                Vec4 vec0;
                Vec4 vec1;
                Vec4 vec2;
                try {
                    vec0 = vertices_transformed.at(triangle.vertexIds[0]);
                    vec1 = vertices_transformed.at(triangle.vertexIds[1]);
                    vec2 = vertices_transformed.at(triangle.vertexIds[2]);
                }
                catch (out_of_range &ex) {
                    continue;
                }

                int x0, y0, x1, y1, x2, y2;
                Color c0, c1, c2;
                x0 = vec0.x;
                y0 = vec0.y;
                c0 = *(this->colorsOfVertices[vec0.colorId - 1]);
                x1 = vec1.x;
                y1 = vec1.y;
                c1 = *(this->colorsOfVertices[vec1.colorId - 1]);
                x2 = vec2.x;
                y2 = vec2.y;
                c2 = *(this->colorsOfVertices[vec2.colorId - 1]);


                line_eq f01(x0, x1, y0, y1);
                line_eq f12(x1, x2, y1, y2);
                line_eq f20(x2, x0, y2, y0);
                int xmin = min(x0, min(x1, x2));
                int ymin = min(y0, min(y1, y2));
                int xmax = max(x0, max(x1, x2));
                int ymax = max(y0, max(y1, y2));
                double alpha1 = 1.0 / (x0 * (f12.y0 - f12.y1) + y0 * (f12.x1 - f12.x0) + f12.x0 * f12.y1 - f12.y0 * f12.x1);
                double beta1 = 1.0 / (x1 * (f20.y0 - f20.y1) + y1 * (f20.x1 - f20.x0) + f20.x0* f20.y1 - f20.y0 * f20.x1);
                double gamma1 = 1.0 / (x2 * (f01.y0 - f01.y1) + y2 * (f01.x1 - f01.x0) + f01.x0 * f01.y1 - f01.y0 * f01.x1);

                for (int y = ymin; y <= ymax; y++) {
                    for (int x = xmin; x <= xmax; x++) {
                        if (x < 0 || y < 0 ||
                            x >= camera->horRes || y >= camera->verRes) {
                            continue;
                        }
                        double alpha = (x * (f12.y0 - f12.y1) + y * (f12.x1 - f12.x0) + f12.x0 * f12.y1 - f12.y0 * f12.x1) * alpha1;
                        double beta = (x * (f20.y0 - f20.y1) + y * (f20.x1 - f20.x0) + f20.x0* f20.y1 - f20.y0 * f20.x1)* beta1;
                        double gamma = (x * (f01.y0 - f01.y1) + y * (f01.x1 - f01.x0) + f01.x0 * f01.y1 - f01.y0 * f01.x1) * gamma1;
                        // Barycentric coordinates
                        /*double alpha2 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) /
                                       ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
                        double beta2 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) /
                                      ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
                        double gamma2 = 1.0 - alpha2 - beta2;*/
                        if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                        //if (isPointInsideTriangle(x+0.5f, y+0.5f, x0,y0,x1,y1,x2,y2)){
                            // Barycentric coordinates
                           /* double alpha2 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) /
                                           ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
                            double beta2 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) /
                                          ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
                            double gamma2 = 1.0 - alpha2 - beta2;*/
                            double interpolatedDepth = alpha * vec0.z + beta * vec1.z + gamma * vec2.z;
                            if (interpolatedDepth <= depthBuffer[y][x]){
                            Color c = c0 * alpha + c1 * beta + c2 * gamma;

                            image[x][y] = c.cround();
                            depthBuffer[y][x] = interpolatedDepth;

                        }


                        }}
                        /*if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                           // double interpolatedDepth = alpha * vec0.z + beta * vec1.z + gamma * vec2.z;
                           // if (interpolatedDepth <= depthBuffer[y][x]){
                                Color c = c0 * alpha + c1 * beta + c2 * gamma;

                                image[x][y] = c.cround();
                               // depthBuffer[y][x] = interpolatedDepth;

                            }*/
                            /*Color c = c0 * alpha + c1 * beta + c2 * gamma;

                            image[x][y] = c.cround();*/
                        }
                    }
                }
            }
        }



