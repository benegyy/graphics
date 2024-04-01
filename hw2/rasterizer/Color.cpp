#include <iomanip>
#include "Color.h"
#include <cmath>
#include <iostream>

Color::Color() {
    this->r = 0;
    this->g = 0;
    this->b = 0;
}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

std::ostream &operator<<(std::ostream &os, const Color &c)
{
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
Color Color::operator-(Color c)
{
    Color c1;
    c1.r = this->r - c.r;
    c1.g = this->g - c.g;
    c1.b = this->b - c.b;

    return c1;
}

Color Color::operator/(int n)
{
    Color c1;
    c1.r = this->r / n;
    c1.g = this->g / n;
    c1.b = this->b / n;

    return c1;
}


Color Color::operator+=(Color c)
{
    this->r = this->r + c.r;
    this->g = this->g + c.g;
    this->b = this->b + c.b;

    return *this;
}
Color Color::operator+(Color c)
{
    Color c1;
    c1.r = this->r + c.r;
    c1.g = this->g + c.g;
    c1.b = this->b + c.b;

    return c1;
}

Color Color::cround()
{
    Color c1;
    c1.r = round(this->r);
    c1.g = round(this->g);
    c1.b = round(this->b);
    return c1;
}

Color Color::operator*(double x){
    Color c1;
    c1.r = this->r * x;
    c1.g = this->g * x;
    c1.b = this->b * x;

    return c1;
}