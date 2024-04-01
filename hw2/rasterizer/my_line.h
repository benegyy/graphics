//
// Created by Bengisu on 12/14/2023.
//

#ifndef CODE_TEMPLATE_MY_LINE_H
#define CODE_TEMPLATE_MY_LINE_H
#include "Vec4.h"
#include "Scene.h"
#include "Mesh.h"
#include "Matrix4.h"
#include "Helpers.h"
class line_eq
{
private:
public:
    int x0,x1,y0,y1;
    line_eq();
    line_eq(int, int, int, int);
    int get_line(int, int);
};
class my_line{
public:
    Vec4 p1,p2;
    Vec4 d;

    my_line();
    my_line(Vec4 p1, Vec4 p2);
    void setD();
};
void line_algorithm(const my_line& line, vector< vector<Color> >& image, const vector< Color* >& colors);

#endif //CODE_TEMPLATE_MY_LINE_H
