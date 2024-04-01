#ifndef __VEC4_H__
#define __VEC4_H__
#define NO_COLOR -1
#include <ostream>
class Vec4
{
public:
    double x, y, z, t;
    int colorId;

    Vec4();
    Vec4(double x, double y, double z, double t);
    Vec4(double x, double y, double z, double t, int colorId);
    Vec4(const Vec4 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec4 &v);
    Vec4 operator=(const Vec4& other);
    Vec4 operator-(const Vec4 p1);
    Vec4 operator+(const Vec4 p1);

};

#endif