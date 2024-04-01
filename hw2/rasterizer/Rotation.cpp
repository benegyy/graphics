#include <iomanip>
#include "Rotation.h"
#include <iostream>
#include "Vec3.h"
#include "Helpers.h"
#include <cmath>
Rotation::Rotation() {
    this->rotationId = -1;
    this->angle = 0;
    this->ux = 0;
    this->uy = 0;
    this->uz = 0;
}

Rotation::Rotation(int rotationId, double angle, double x, double y, double z)
{
    this->rotationId = rotationId;
    this->angle = angle;
    this->ux = x;
    this->uy = y;
    this->uz = z;
}

std::ostream &operator<<(std::ostream &os, const Rotation &r)
{
    os << std::fixed << std::setprecision(3) << "Rotation " << r.rotationId << " => [angle=" << r.angle << ", " << r.ux << ", " << r.uy << ", " << r.uz << "]";
    return os;
}
Matrix4 Rotation::rotate_matrix()
{
    Matrix4 result;
    Vec3 u;
    Vec3 v;

    u.x = ux;
    u.y = uy;
    u.z = uz;

    double x = normalizeVec3(u).x;
    double y = normalizeVec3(u).y;
    double z = normalizeVec3(u).z;

    u.x = x;
    u.y = y;
    u.z = z;

    double min_element = u.x;
    char element_index = 'x';

    if(abs(u.y) < abs(min_element)){
        min_element = u.y;
        element_index = 'y';
    }
    if(abs(u.z) < abs(min_element)){
        min_element = u.z;
        element_index = 'z';
    }

    if(element_index == 'x'){
        v.x = 0.0;
        v.y = (-1.0)*u.z;
        v.z = u.y;
    }
    else if(element_index == 'y'){
        v.x = (-1.0)*u.z;
        v.y = 0.0;
        v.z = u.x;
    }
    else{
        v.x = (-1.0)*u.y;
        v.y = u.x;
        v.z = 0.0;
    }


    x = normalizeVec3(v).x;
    y = normalizeVec3(v).y;
    z = normalizeVec3(v).z;

    v.x = x;
    v.y = y;
    v.z = z;
    Vec3 w= crossProductVec3(u, v);
    Matrix4 M;
    Matrix4 M_inverse;

    M_inverse.values[0][0] = u.x;
    M_inverse.values[0][1] = v.x;
    M_inverse.values[0][2] = w.x;
    M_inverse.values[0][3] = 0.0;

    M_inverse.values[1][0] = u.y;
    M_inverse.values[1][1] = v.y;
    M_inverse.values[1][2] = w.y;
    M_inverse.values[1][3] = 0.0;

    M_inverse.values[2][0] = u.z;
    M_inverse.values[2][1] = v.z;
    M_inverse.values[2][2] = w.z;
    M_inverse.values[2][3] = 0.0;

    M_inverse.values[3][0] = 0.0;
    M_inverse.values[3][1] = 0.0;
    M_inverse.values[3][2] = 0.0;
    M_inverse.values[3][3] = 1.0;



    M.values[0][0] = u.x;
    M.values[0][1] = u.y;
    M.values[0][2] = u.z;
    M.values[0][3] = 0.0;

    M.values[1][0] = v.x;
    M.values[1][1] = v.y;
    M.values[1][2] = v.z;
    M.values[1][3] = 0.0;

    M.values[2][0] = w.x;
    M.values[2][1] = w.y;
    M.values[2][2] = w.z;
    M.values[2][3] = 0.0;

    M.values[3][0] = 0.0;
    M.values[3][1] = 0.0;
    M.values[3][2] = 0.0;
    M.values[3][3] = 1.0;


    double theta = (angle*M_PI)/180.0;
    Matrix4 Rx_theta;

    Rx_theta.values[0][0] = 1.0;
    Rx_theta.values[0][1] = 0.0;
    Rx_theta.values[0][2] = 0.0;
    Rx_theta.values[0][3] = 0.0;

    Rx_theta.values[1][0] = 0.0;
    Rx_theta.values[1][1] = cos(theta);
    Rx_theta.values[1][2] = (-1.0)*sin(theta);
    Rx_theta.values[1][3] = 0.0;

    Rx_theta.values[2][0] = 0.0;
    Rx_theta.values[2][1] = sin(theta);
    Rx_theta.values[2][2] = cos(theta);
    Rx_theta.values[2][3] = 0.0;

    Rx_theta.values[3][0] = 0.0;
    Rx_theta.values[3][1] = 0.0;
    Rx_theta.values[3][2] = 0.0;
    Rx_theta.values[3][3] = 1.0;


    Matrix4 Rx_times_M;

    Rx_times_M=multiplyMatrixWithMatrix( Rx_theta, M);
    result=multiplyMatrixWithMatrix( M_inverse, Rx_times_M);
    return result;
}