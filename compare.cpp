#include <iostream>
#include <random>

#define _USE_MATH_DEFINES
#include <cmath>

#include <cassert>

#include "geo.h"
#include "bvhbox.h"
#include "bvhsphere.h"
#include "kdtree.h"
#include "quadtree.h"
#include "vptree.h"
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

void initbunny(int dx, int dy, std::vector<geo::Point>& points)
{
    points.clear();
    std::ifstream file("bunny.dat");
    geo::Point p;
    while(file)
    {
        file >> p.x >> p.y;
        points.push_back(geo::Point(p.x*dx, p.y*dy));
    }
}

void initdiag(int dx, int dy, std::vector<geo::Point>& points)
{
    int nb = 126;
    points.clear();
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        float x = rand1D() * dx;
        points[i] = geo::Point(x,x);
    }
}

static float gaussian(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;
    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

void initgauss(int dx, int dy, std::vector<geo::Point>& points)
{
    int nb = 126;
    points.clear();
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        float x = rand1D();
        float y = gaussian(3.f*x, 1.5f, 0.4f);
        points[i] = geo::Point(x*dx,y*dy);
    }
}

void initsquare(int dx, int dy, std::vector<geo::Point>& points)
{
    int nb = 126;
    points.clear();
    points.resize(nb);
    for(int i = 0; i < nb/2; ++i)
    {
        float x = rand1D();
        points[2*i] = geo::Point(x*dx,x*x*dy);
        points[2*i+1] = geo::Point(x*dx,(1.f-x*x)*dy);
    }
}

int main(void)
{   
    int dx = 1000, dy = 1000;
    svg::Dimensions dimensions(dx, dy);
    svg::Layout layout(dimensions, svg::Layout::BottomLeft);
    svg::Polygon bg(svg::Fill(svg::Color(255, 255, 255)), svg::Stroke());
    bg << svg::Point(0, 0) << svg::Point(dimensions.width, 0) << svg::Point(dimensions.width, dimensions.height) << svg::Point(0, dimensions.height);

    // init pointset
    std::vector<geo::Point> v;
    initbunny(dx,dy,v);
    // initdiag(dx,dy,v);
    // initgauss(dx,dy,v);
    // initsquare(dx,dy,v);

    // bvh sphere tree
    svg::Document d1("bvhsphere.svg", layout);
    d1 << bg;
    bvhsphere::Tree bsh(v);
    bvhsphere::draw(bsh, 0, bsh.root, d1);
    d1.save();

    // bvh box tree
    svg::Document d2("bvhbox.svg", layout);
    d2 << bg;
    bvhbox::Tree bbh(v);
    bvhbox::draw(bbh, 0, bbh.root, d2);
    d2.save();

    // quad tree
    svg::Document d3("quadtree.svg", layout);
    d3 << bg;
    quadtree::Tree qt(v, geo::Box(geo::Point(0,0), geo::Point(dx, dy)));
    quadtree::draw(qt, 0, qt.root, d3);
    d3.save();

    // kd tree
    svg::Document d4("kdtree.svg", layout);
    d4 << bg;
    kdtree::Tree kd(v);
    kdtree::draw(kd, 0, kd.root, geo::Point(0,0), geo::Point(dx, dy), d4);
    d4.save();
    
    // vp tree 
    svg::Document d5("vptree.svg", layout);
    d5 << bg;
    vptree::Tree vp(v);
    // init canvas cell
    cavc::Polyline<double> bounds;
    bounds.addVertex(1, 1, 0); 
    bounds.addVertex(dx-1, 1, 0); 
    bounds.addVertex(dx-1, dy-1, 0); 
    bounds.addVertex(1, dy-1, 0);
    bounds.isClosed() = true;
    vptree::Cell cell;
    cell.shape.push_back(bounds);
    // init edge structure
    std::vector<svg::CavcPoly::Edge> buffer;
    vptree::draw(vp, 0, vp.root, cell, buffer, d5);
    d5.save();

    return 0;
}
