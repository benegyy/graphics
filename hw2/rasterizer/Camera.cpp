#include <iomanip>
#include "Camera.h"
#include <string>
#include <iostream>
#include <iomanip>
#include "Matrix4.h"
#include "Helpers.h"


Camera::Camera() {}

Camera::Camera(int cameraId,
               int projectionType,
               Vec3 position, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               std::string outputFilename)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->position = position;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFilename = outputFilename;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->position = other.position;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFilename = other.outputFilename;
}

std::ostream &operator<<(std::ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << std::fixed << std::setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.position << " gaze: " << c.gaze << std::endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << std::endl
       << std::fixed << std::setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << std::endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFilename;

    return os;
}

Matrix4 Camera::camera_transfromation()
{

    Matrix4 result;
    result.values[0][0] = u.x;
    result.values[0][1] = u.y;
    result.values[0][2] = u.z;
    result.values[0][3] = -dotProductVec3(u, position);
    result.values[1][0] = v.x;
    result.values[1][1] = v.y;
    result.values[1][2] = v.z;
    result.values[1][3] = -dotProductVec3(v, position);
    result.values[2][0] = w.x;
    result.values[2][1] = w.y;
    result.values[2][2] = w.z;
    result.values[2][3] = -dotProductVec3(w, position);
    result.values[3][0] = 0.0;
    result.values[3][1] = 0.0;
    result.values[3][2] = 0.0;
    result.values[3][3] = 1;

    return result;
}

Matrix4 Camera::perspective_and_orth( )
{        Matrix4 resultf;

    if(!projectionType){
        Matrix4 result;
    result.values[0][0] = 2 / (right-left);
    result.values[0][1] = 0.0;
    result.values[0][2] = 0.0;
    result.values[0][3] = -(right+left) / (right-left);
    result.values[1][0] = 0;
    result.values[1][1] = 2 / (top-bottom);
    result.values[1][2] = 0;
    result.values[1][3] =-(top+bottom) / (top-bottom);
    result.values[2][0] =0;
    result.values[2][1] = 0;
    result.values[2][2] = -2 / (far-near);
    result.values[2][3] = -(far+near) / (far-near);
    result.values[3][0] =0;
    result.values[3][1] = 0;
    result.values[3][2] = 0;
    result.values[3][3] = 1;
    return result;
    }

    else{
        Matrix4 result2;
        result2.values[0][0] = (2.0)*near/(right - left);
        result2.values[0][1] = 0.0;
        result2.values[0][2] = (right + left)/(right - left);
        result2.values[0][3] = 0.0;

        result2.values[1][0] = 0.0;
        result2.values[1][1] = (2.0)*near / (top - bottom);
        result2.values[1][2] = (top + bottom)/(top - bottom);
        result2.values[1][3] = 0.0;

        result2.values[2][0] = 0.0;
        result2.values[2][1] = 0.0;
        result2.values[2][2] = (-1.0)*(far + near) / (far - near);
        result2.values[2][3] = (-2.0)*far*near / (far - near);

        result2.values[3][0] = 0.0;
        result2.values[3][1] = 0.0;
        result2.values[3][2] = -1.0;
        result2.values[3][3] = 0.0;
        return result2;
    }
}

Matrix4 Camera::computeMperMcam()
{
    return multiplyMatrixWithMatrix(perspective_and_orth(),
                                    camera_transfromation());
}

Matrix4 Camera::viewport_transformation()
{
    Matrix4 result;

    result.values[0][0] = horRes/2.0;
    result.values[0][1] = 0.0;
    result.values[0][2] = 0.0;
    result.values[0][3] = (horRes-1.0)/2.0;

    result.values[1][0] = 0.0;
    result.values[1][1] = verRes/2.0;
    result.values[1][2] = 0.0;
    result.values[1][3] = (verRes-1.0)/2.0;

    result.values[2][0] = 0.0;
    result.values[2][1] = 0.0;
    result.values[2][2] = 0.5;
    result.values[2][3] = 0.5;

    result.values[3][0] = 0;
    result.values[3][1] = 0;
    result.values[3][2] = 0;
    result.values[3][3] = 1;

    return result;
}