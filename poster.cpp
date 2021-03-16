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
    svg::Dimensions dimensions(4000, 5000);
    svg::Layout layout(dimensions, svg::Layout::BottomLeft);
    svg::Polygon bg(svg::Fill(svg::Color(255, 255, 255)), svg::Stroke());
    bg << svg::Point(0, 0) << svg::Point(dimensions.width, 0) << svg::Point(dimensions.width, dimensions.height) << svg::Point(0, dimensions.height);
    svg::Document doc("poster.svg", layout);
    doc << bg;

    if( 1 )
    {
        int sx = 800, sy = 800;
        geo::Point offset;
        std::vector<geo::Point> v;
        initbunny(sx, sy, v);
        std::vector<geo::Point> vcpy(v.size());
        // bunny bvh sphere
        offset = geo::Point(100,100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhsphere::Tree t1(vcpy);
        bvhsphere::draw(t1, 0, t1.root, doc);
        // bunny bvh box
        offset = geo::Point(100,1100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhbox::Tree t2(vcpy);
        bvhbox::draw(t2, 0, t2.root, doc);
        // bunny kd tree
        offset = geo::Point(100,2100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        kdtree::Tree t3(vcpy);
        kdtree::draw(t3, 0, t3.root, geo::Point(100,2100), geo::Point(900, 2900), doc);
        // bunny quad tree
        offset = geo::Point(100,3100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        quadtree::Tree t4(vcpy, geo::Box(geo::Point(100,3100), geo::Point(900,3900)));
        quadtree::draw(t4, 0, t4.root, doc);
        // bunny vp tree
        offset = geo::Point(100,4100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        vptree::Tree t5(vcpy);
        cavc::Polyline<double> bounds;
        bounds.addVertex(100, 4100, 0); 
        bounds.addVertex(100, 4900, 0); 
        bounds.addVertex(900, 4900, 0); 
        bounds.addVertex(900, 4100, 0);
        bounds.isClosed() = true;
        vptree::Cell cell;
        cell.shape.push_back(bounds);
        std::vector<svg::CavcPoly::Edge> buffer;
        vptree::draw(t5, 0, t5.root, cell, buffer, doc);
    }
    
    if( 1 )
    {
        int sx = 800, sy = 800;
        geo::Point offset;
        std::vector<geo::Point> v;
        initdiag(sx, sy, v);
        std::vector<geo::Point> vcpy(v.size());
        // line bvh sphere
        offset = geo::Point(1100,100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhsphere::Tree t1(vcpy);
        bvhsphere::draw(t1, 0, t1.root, doc);
        // line bvh box
        offset = geo::Point(1100,1100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhbox::Tree t2(vcpy);
        bvhbox::draw(t2, 0, t2.root, doc);
        // line kd tree
        offset = geo::Point(1100,2100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        kdtree::Tree t3(vcpy);
        kdtree::draw(t3, 0, t3.root, geo::Point(1100,2100), geo::Point(1900, 2900), doc);
        // line quad tree
        offset = geo::Point(1100,3100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        quadtree::Tree t4(vcpy, geo::Box(geo::Point(1100,3100), geo::Point(1900,3900)));
        quadtree::draw(t4, 0, t4.root, doc);
        // line vp tree
        offset = geo::Point(1100,4100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        vptree::Tree t5(vcpy);
        cavc::Polyline<double> bounds;
        bounds.addVertex(1100, 4100, 0); 
        bounds.addVertex(1100, 4900, 0); 
        bounds.addVertex(1900, 4900, 0); 
        bounds.addVertex(1900, 4100, 0);
        bounds.isClosed() = true;
        vptree::Cell cell;
        cell.shape.push_back(bounds);
        std::vector<svg::CavcPoly::Edge> buffer;
        vptree::draw(t5, 0, t5.root, cell, buffer, doc);
    }
    
    if( 1 )
    {
        int sx = 800, sy = 800;
        geo::Point offset;
        std::vector<geo::Point> v;
        initsquare(sx, sy, v);
        std::vector<geo::Point> vcpy(v.size());
        // squares bvh sphere
        offset = geo::Point(2100,100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhsphere::Tree t1(vcpy);
        bvhsphere::draw(t1, 0, t1.root, doc);
        // squares bvh box
        offset = geo::Point(2100,1100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhbox::Tree t2(vcpy);
        bvhbox::draw(t2, 0, t2.root, doc);
        // squares kd tree
        offset = geo::Point(2100,2100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        kdtree::Tree t3(vcpy);
        kdtree::draw(t3, 0, t3.root, geo::Point(2100,2100), geo::Point(2900, 2900), doc);
        // squares quad tree
        offset = geo::Point(2100,3100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        quadtree::Tree t4(vcpy, geo::Box(geo::Point(2100,3100), geo::Point(2900,3900)));
        quadtree::draw(t4, 0, t4.root, doc);
        // squares vp tree
        offset = geo::Point(2100,4100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        vptree::Tree t5(vcpy);
        cavc::Polyline<double> bounds;
        bounds.addVertex(2100, 4100, 0); 
        bounds.addVertex(2100, 4900, 0); 
        bounds.addVertex(2900, 4900, 0); 
        bounds.addVertex(2900, 4100, 0);
        bounds.isClosed() = true;
        vptree::Cell cell;
        cell.shape.push_back(bounds);
        std::vector<svg::CavcPoly::Edge> buffer;
        vptree::draw(t5, 0, t5.root, cell, buffer, doc);
    }
    
    if( 1 )
    {
        int sx = 800, sy = 800;
        geo::Point offset;
        std::vector<geo::Point> v;
        initgauss(sx, sy, v);
        std::vector<geo::Point> vcpy(v.size());
        // gaussian bvh sphere
        offset = geo::Point(3100,100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhsphere::Tree t1(vcpy);
        bvhsphere::draw(t1, 0, t1.root, doc);
        // gaussian bvh box
        offset = geo::Point(3100,1100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        bvhbox::Tree t2(vcpy);
        bvhbox::draw(t2, 0, t2.root, doc);
        // gaussian kd tree
        offset = geo::Point(3100,2100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        kdtree::Tree t3(vcpy);
        kdtree::draw(t3, 0, t3.root, geo::Point(3100,2100), geo::Point(3900, 2900), doc);
        // gaussian quad tree
        offset = geo::Point(3100,3100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        quadtree::Tree t4(vcpy, geo::Box(geo::Point(3100,3100), geo::Point(3900,3900)));
        quadtree::draw(t4, 0, t4.root, doc);
        // gaussian vp tree
        offset = geo::Point(3100,4100);
        for(int i = 0; i < (int)v.size(); ++i)
            vcpy[i] = v[i]+offset;
        vptree::Tree t5(vcpy);
        cavc::Polyline<double> bounds;
        bounds.addVertex(3100, 4100, 0); 
        bounds.addVertex(3100, 4900, 0); 
        bounds.addVertex(3900, 4900, 0); 
        bounds.addVertex(3900, 4100, 0);
        bounds.isClosed() = true;
        vptree::Cell cell;
        cell.shape.push_back(bounds);
        std::vector<svg::CavcPoly::Edge> buffer;
        vptree::draw(t5, 0, t5.root, cell, buffer, doc);
    }
    doc.save();


    return 0;
}
