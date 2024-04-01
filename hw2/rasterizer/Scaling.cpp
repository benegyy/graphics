#include <iomanip>
#include "Scaling.h"
#include <iostream>
#include "Matrix4.h"
Scaling::Scaling() {
    this->scalingId = -1;
    this->sx = 0;
    this->sy = 0;
    this->sz = 0;
}

Scaling::Scaling(int scalingId, double sx, double sy, double sz)
{
    this->scalingId = scalingId;
    this->sx = sx;
    this->sy = sy;
    this->sz = sz;
}

std::ostream &operator<<(std::ostream &os, const Scaling &s)
{
    os << std::fixed << std::setprecision(3) << "Scaling " << s.scalingId << " => [" << s.sx << ", " << s.sy << ", " << s.sz << "]";

    return os;
}
Matrix4 Scaling::scale_matrix()
{
    Matrix4 result;
    result.values[0][0] = sx;
    result.values[0][1] = 0;
    result.values[0][2] = 0;
    result.values[0][3] = 0;

    result.values[1][0] = 0;
    result.values[1][1] = sy;
    result.values[1][2] = 0;
    result.values[1][3] = 0;

    result.values[2][0] = 0;
    result.values[2][1] = 0;
    result.values[2][2] = sz;
    result.values[2][3] = 0;

    result.values[3][0] = 0;
    result.values[3][1] = 0;
    result.values[3][2] = 0;
    result.values[3][3] = 1.0;

    return  result;
}