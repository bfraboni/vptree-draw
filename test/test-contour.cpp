#include <iostream>
#include <random>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cassert>
#include <stack>

#include "geo.h"
#include "bvhbox.h"
#include "bvhsphere.h"
#include "kdtree.h"
#include "quadtree.h"
#include "vptree.h"
#include "bvptree.h"
#include "draw.h"

static float rand1D() 
{   
    using G = std::default_random_engine;    
    using D = std::uniform_real_distribution<float>; 
    static const std::size_t seed = 123456789;   
    static G gen(seed);
    static D dist(0.f, 1.f);
    return dist(gen);
}


static float bregmanKL(const geo::Point& a, const geo::Point & b)
{
    return a.x*std::log(a.x/b.x)+a.y*std::log(a.y/b.y)+b.x-a.x+b.y-a.y;
}

// A utility function to find next to top in a stack
geo::Point nextToTop(std::stack<geo::Point> &S)
{
    geo::Point p = S.top();
    S.pop();
    geo::Point res = S.top();
    S.push(p);
    return res;
}

int main(int argc, char * argv[])
{   
    int dx = 1000, dy = 1000;
    svg::Dimensions dimensions(dx, dy);
    svg::Layout layout(dimensions, svg::Layout::BottomLeft);
    svg::Polygon bg(svg::Fill(svg::Color(255, 255, 255)), svg::Stroke());
    bg << svg::Point(0, 0) << svg::Point(dimensions.width, 0) << svg::Point(dimensions.width, dimensions.height) << svg::Point(0, dimensions.height);
    svg::Document doc("contour.svg", layout);
    doc << bg;

    geo::Point p(dx/2,dy/2);
    float tau = argc > 1 ? std::atof(argv[1]) : 1 ;

    // cf https://stackoverflow.com/questions/4313992/methods-for-implementing-contour-plotting
    std::vector<geo::Point> v;
    for(int i = -400; i < 400; i+=4)
    for(int j = -400; j < 400; j+=4)
    {
        {
            geo::Point a = p + geo::Point(i,j);
            geo::Point b = p + geo::Point(i+2,j);
            float fa = bregmanKL(a,p);
            float fb = bregmanKL(b,p);

            if( (fa-tau) * (fb-tau) < 0 )
            {
                float t = fa-tau > 0 ? (fa-tau) / (fa-fb) : (fb-tau) / (fb-fa);
                v.push_back(lerp(a,b,t));
                doc << svg::Circle(svg::Point(v.back().x, v.back().y), 4, svg::Fill(svg::Color(0,0,0)), svg::Stroke());
            }
        }

        {
            geo::Point a = p + geo::Point(i,j);
            geo::Point b = p + geo::Point(i,j+2);
            float fa = bregmanKL(a,p);
            float fb = bregmanKL(b,p);

            if( (fa-tau) * (fb-tau) < 0 )
            {
                float t = fa-tau > 0 ? (fa-tau) / (fa-fb) : (fb-tau) / (fb-fa);
                v.push_back(lerp(a,b,t));
                doc << svg::Circle(svg::Point(v.back().x, v.back().y), 4, svg::Fill(svg::Color(0,0,0)), svg::Stroke());
            }
        }
    }

    // cf : https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/

    // Find the bottommost point
    int ymin = v[0].y, min = 0;
    for (int i = 1; i < v.size(); i++)
    {
        if((v[i].y < ymin) || (ymin == v[i].y && v[i].x < v[min].x))
            ymin = v[i].y, min = i;
    }

    // Place the bottom-most point at first position
    std::swap(v[0], v[min]);

    // Sort n-1 points with respect to the first point.
    // A point p1 comes before p2 in sorted output if p2
    // has larger polar angle (in counterclockwise
    // direction) than p1
    geo::Point p0 = v[0];
    auto compare = [&p0](const geo::Point& a, const geo::Point& b) 
    {
        float val = (a.y - p0.y) * (b.x - a.x) - (a.x - p0.x) * (b.y - a.y);
        if( val < 0 ) return true;
        if( val > 0 ) return false;
        return (distance2(p0, b) >= distance2(p0, a));
    };

    std::sort(v.begin(), v.end(), compare);

    // Create an empty stack and push first three points
    // to it.
    std::stack<geo::Point> stack;
    stack.push(v[0]);
    stack.push(v[1]);
    stack.push(v[2]);

    auto orientation = [](const geo::Point& p, const geo::Point& a, const geo::Point& b)
    {
        return (a.y - p.y) * (b.x - a.x) - (a.x - p.x) * (b.y - a.y); // colinear = 0, clock > 0 or counterclock wise < 0
    };

    // Process remaining n-3 points
    for (int i = 3; i < v.size(); i++)
    {
        // Keep removing top while the angle formed by
        // points next-to-top, top, and points[i] makes
        // a non-left turn
        while(stack.size() > 1 && orientation(nextToTop(stack), stack.top(), v[i]) >= 0){
            stack.pop();
        }
        stack.push(v[i]);
    }

    // Now stack has the output points, print contents of stack
    svg::Polygon hull(svg::Fill(), svg::Stroke(2, svg::Color(255, 0, 0)));
    while (!stack.empty())
    {
        geo::Point p = stack.top();
        stack.pop();
        hull << svg::Point(p.x, p.y);
    }
    doc << hull;
    doc.save();
    return 0;
}
