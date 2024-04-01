#include "my_line.h"
#include <cmath>
#include <vector>

my_line::my_line(Vec4 point1, Vec4 point2)
{
    this->p1 = point1;
    this->p2 = point2;

}

void my_line::setD(){
    this->d =  p2-p1;
}
line_eq::line_eq(int x0, int x1, int y0, int y1)
{
    this->x0 = x0;
    this->y0 = y0;
    this->x1 = x1;
    this->y1 = y1;
}

int line_eq::get_line(int x, int y){
    return x*(y0-y1)+y*(x1-x0)+x0*y1-y0*x1;
}


void line_algorithm(const my_line& line, vector< vector<Color> >& image, const vector< Color* >& colors)
{
    double slope;
    int x0, x1, y0, y1;
    Color color0, color1;
    Color color(0, 0, 0);

    slope = (line.p1.y - line.p2.y)/(line.p1.x - line.p2.x);
    if(0 <= slope)
    {
        if(slope <1)
        {
            if(line.p1.x < line.p2.x)
            {
                x0 = round(line.p1.x);
                y0 = round(line.p1.y);
                color0 = *colors[line.p1.colorId-1];

                x1 = round(line.p2.x);
                y1 = round(line.p2.y);
                color1 = *colors[line.p2.colorId-1];
            }
            else
            {
                x0 = round(line.p2.x);
                y0 = round(line.p2.y);
                color0 = *colors[line.p2.colorId-1];

                x1 = round(line.p1.x);
                y1 = round(line.p1.y);
                color1 = *colors[line.p1.colorId-1];
            }

            int y = y0;
            int d = 2*(y0 - y1) + (x1 - x0);
            Color color_a = color0;
            Color color_b = (color1 - color0) / (x1 - x0);
            for(int x = x0; x <= x1; x++)
            {
                image[x][y] = color_a.cround();

                if(d < 0){
                    y++;
                    d += 2*((y0-y1)+(x1-x0));
                }
                else{
                    d += 2*(y0-y1);
                }
                color_a += color_b;
            }
        }
        else
        {
            if(line.p1.x < line.p2.x)//
            {
                x0 = round(line.p1.x);
                y0 = round(line.p1.y);
                color0 = *colors[line.p1.colorId-1];

                x1 = round(line.p2.x);
                y1 = round(line.p2.y);
                color1 = *colors[line.p2.colorId-1];
            }
            else
            {
                x0 = round(line.p2.x);
                y0 = round(line.p2.y);
                color0 = *colors[line.p2.colorId-1];

                x1 = round(line.p1.x);
                y1 = round(line.p1.y);
                color1 = *colors[line.p1.colorId-1];
            }

            int x = x0;
            int d = (y0 - y1) + 2*(x1 - x0);
            Color c = color0;
            Color dc = (color1 - color0)/(y1 - y0);
            for(int y = y0; y <= y1; y++)
            {
                image[x][y] = c.cround();
                if(d < 0)
                {
                    d += 2*(x1 - x0);
                }
                else
                {
                    x++;
                    d += 2*((y0 - y1) + (x1 - x0 ));
                }
                c += dc;
            }
        }
    }
    else
    {
        if( -1 >= slope)
        {
            if(line.p1.y < line.p2.y)
            {
                x0 = round(line.p1.x);
                y0 = round(line.p1.y);
                color0 = *colors[line.p1.colorId-1];

                x1 = round(line.p2.x);
                y1 = round(line.p2.y);
                color1 = *colors[line.p2.colorId-1];
            }
            else
            {
                x0 = round(line.p2.x);
                y0 = round(line.p2.y);
                color0 = *colors[line.p2.colorId-1];

                x1 = round(line.p1.x);
                y1 = round(line.p1.y);
                color1 = *colors[line.p1.colorId-1];
            }

            int x = x0;
            int d = -(y0 - y1) + 2*(x1 - x0);
            Color c = color0;
            Color dc = (color1 - color0)/(y1 - y0);

            for(int y = y0; y<= y1; y++)
            {
                image[x][y] = c.cround();

                if(d < 0){
                    x--;
                    d += 2*(-(y0 - y1) + (x1 - x0) );
                }
                else
                {
                    d += 2*(x1 - x0);
                }
                c += dc;
            }
        }
        else
        {

            if(line.p1.y < line.p2.y)
            {
                x0 = round(line.p1.x);
                y0 = round(line.p1.y);
                color0 = *colors[line.p1.colorId-1];

                x1 = round(line.p2.x);
                y1 = round(line.p2.y);
                color1 = *colors[line.p2.colorId-1];
            }
            else
            {
                x0 = round(line.p2.x);
                y0 = round(line.p2.y);
                color0 = *colors[line.p2.colorId-1];

                x1 = round(line.p1.x);
                y1 = round(line.p1.y);
                color1 = *colors[line.p1.colorId-1];
            }


            int y = y0;
            int d = -2*(y0 - y1) + (x1 - x0);
            Color c = color0;
            Color dc = (color1 - color0)/(x0-x1);

            for(int x = x0; x >= x1; x--)
            {
                image[x][y] = c.cround();

                if( d < 0){
                    d += -2 *(y0 - y1);
                }
                else{
                    y++;
                    d += 2* (-(y0 - y1) + (x1 -x0));
                }
                c += dc;
            }
        }
    }
}
