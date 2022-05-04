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
#include "bvptree.h"
#include "draw.h"

int main(int argc, char * argv[])
{   
    int dx = 1000, dy = 1000;
    svg::Dimensions dimensions(dx, dy);
    svg::Layout layout(dimensions, svg::Layout::BottomLeft);
    svg::Polygon bg(svg::Fill(svg::Color(255, 255, 255)), svg::Stroke());
    bg << svg::Point(0, 0) << svg::Point(dimensions.width, 0) << svg::Point(dimensions.width, dimensions.height) << svg::Point(0, dimensions.height);

    // bregman ball size
    geo::Point c(200, 200);
    float tau = 200.f;
    if( argc > 1 )
    {
        tau = std::atof(argv[3]);
    }

    // old version
    if(0)
    {
        std::vector<geo::Point> vold = geo::ball(c, tau, geo::BregmanKL());
        std::vector<geo::Point> hull = geo::hull(vold);
    
        svg::Polygon poly(svg::Fill(), svg::Stroke(1, svg::Color(0,0,0))); 
        for(auto p : hull)
            poly << svg::Point(p.x, p.y);

        svg::Document doc("klball-old.svg", layout);
    
        doc << bg;
        doc << poly;
        doc.save();
    }

    // new version
    if(0)
    {
        std::vector<geo::Point> vnew = geo::klball(c, tau);

        svg::Polygon poly(svg::Fill(), svg::Stroke(1, svg::Color(0,0,0))); 
        for(auto p : vnew)
            poly << svg::Point(p.x, p.y);

        svg::Document doc("klball-new.svg", layout);
        doc << bg;
        doc << poly;
        doc.save();
    }

    int nbx = 3, nby = 3, nbr = 4;
    // old version
    if(1)
    {
        svg::Document doc("klball-old-demo.svg", layout);
        doc << bg;
        std::vector<geo::Point> vold, hull;
        for(int i = 0; i < nbx; ++i)
        for(int j = 0; j < nby; ++j)
        for(int k = 0; k < nbr; ++k)
        {
            geo::Point c((i+0.5f)*dx/float(nbx),(j+0.5f)*dy/float(nby));
            float tau = (k+1)*25;
            vold = geo::ball(c, tau, geo::BregmanKL());
            hull = geo::hull(vold);

            svg::Color color = rotate(svg::Color(255,120,80), dx*(k+1.f)/10.f);
            svg::Polygon poly(svg::Fill(), svg::Stroke(1, color)); 
            for(auto p : hull)
                poly << svg::Point(p.x, p.y);
        
            doc << poly;
        }
        doc.save();
    }

    // new version
    if(1)
    {
        svg::Document doc("klball-new-demo.svg", layout);
        doc << bg;
        std::vector<geo::Point> vnew;
        for(int i = 0; i < nbx; ++i)
        for(int j = 0; j < nby; ++j)
        for(int k = 0; k < nbr; ++k)
        {
            geo::Point c((i+0.5f)*dx/float(nbx),(j+0.5f)*dy/float(nby));
            float tau = (k+1)*25;
            vnew = geo::klball(c, tau);

            svg::Color color = rotate(svg::Color(255,120,80), dx*(k+1.f)/10.f);
            svg::Polygon poly(svg::Fill(), svg::Stroke(1, color)); 
            for(auto p : vnew)
                poly << svg::Point(p.x, p.y);
        
            doc << poly;
        }
        doc.save();
    }

    return 0;
}
