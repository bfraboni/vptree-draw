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

void initrandom(int dx, int dy, std::vector<geo::Point>& points)
{
    int nb = 100;
    points.clear();
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        float x = rand1D();
        float y = rand1D();
        points[i] = geo::Point(x*dx/2+dx/4, y*dy/2+dy/4);
    }
}

namespace geo
{
    // https://en.wikipedia.org/wiki/Smallest-circle_problem#Linear-time_solutions
    // Sometimes french wikipedia is better :)
    // https://fr.wikipedia.org/wiki/Probl%C3%A8me_du_cercle_minimum
    geo::Sphere minisphere2(std::vector<geo::Point>& points, int i, std::vector<geo::Point> r)
    {

        printf("D %d P ", i);
        for(int j = 0; j < points.size(); ++j){
            if(j==i) printf(" i ");
            printf("(%f, %f) ", points[j].x, points[j].y);
        }
        printf("R ");
        for(int j = 0; j < r.size(); ++j)
            printf("(%f, %f) ", points[j].x, points[j].y);
        printf("\n");

        if(i == points.size() || r.size() == 3)
        {
            // printf("trivial case %d points\n", int(r.size()));
            if( r.size() == 0) return geo::Sphere(geo::Point::zero(), -1);
            if( r.size() == 1) return geo::Sphere(r[0], 0);
            if( r.size() == 2) return geo::Sphere(r[0], r[1]);
            if( r.size() == 3) return geo::Sphere(r[0], r[1], r[2]);
        }
        geo::Sphere s = minisphere2(points, i+1, r);
        if(!s.contains(points[i]))
        {
            r.push_back(points[i]);
            s = minisphere2(points, i+1, r);
        }
        return s;
    }

    geo::Sphere minisphere2(std::vector<geo::Point>& points)
    {
        std::shuffle(points.begin(), points.end(), std::mt19937());
        geo::Sphere s = minisphere2(points, 0, std::vector<geo::Point>());
        assert(s.r >= 0);
        return s;
    }

    // https://en.wikipedia.org/wiki/Smallest-circle_problem#Linear-time_solutions
    // Sometimes french wikipedia is better :)
    // https://fr.wikipedia.org/wiki/Probl%C3%A8me_du_cercle_minimum
    geo::Sphere minisphere3(std::vector<geo::Point>& points, int i, int r)
    {
        // if( i > 6) return geo::Sphere(geo::Point::zero(), -1);
        printf("D %d P ", i);
        for(int j = 0; j < points.size(); ++j){
            if(j==i) printf(" i ");
            printf("(%f, %f) ", points[j].x, points[j].y);
        }
        printf("R ");
        for(int j = 0; j < r; ++j)
            printf("(%f, %f) ", points[j].x, points[j].y);
        printf("\n");

        if(i == points.size() || r+1 == 3)
        {
            // printf("trivial case %d points\n", int(r.size()));
            if( r == -1) return geo::Sphere(geo::Point::zero(), -1);
            if( r == 0) return geo::Sphere(points[0], 0);
            if( r == 1) return geo::Sphere(points[0], points[1]);
            if( r == 2) return geo::Sphere(points[0], points[1], points[2]);
        }
        geo::Sphere s = minisphere3(points, i+1, r);
        if(!s.contains(points[i]))
        {
            s = minisphere3(points, i+1, r+1);
        }
        return s;
    }

    geo::Sphere minisphere3(std::vector<geo::Point>& points)
    {
        std::shuffle(points.begin(), points.end(), std::mt19937());
        geo::Sphere s = minisphere3(points, 0, -1);
        assert(s.r >= 0);
        return s;
    }
}

int main(void)
{   
    int dx = 1000, dy = 1000;
    svg::Dimensions dimensions(dx, dy);
    svg::Layout layout(dimensions, svg::Layout::BottomLeft);
    svg::Polygon bg(svg::Fill(svg::Color(255, 255, 255)), svg::Stroke());
    bg << svg::Point(0, 0) << svg::Point(dimensions.width, 0) << svg::Point(dimensions.width, dimensions.height) << svg::Point(0, dimensions.height);
    svg::Document doc("minisphere.svg", layout);
    doc << bg;


    std::vector<geo::Point> v;
    initrandom(dx,dy,v);
    std::shuffle(v.begin(), v.end(), std::mt19937());
    printf("\nmini2\n");
    geo::Sphere s1 = minisphere2(v);
    printf("\nmini3\n");
    geo::Sphere s = minisphere3(v);

    printf("\nmini2 (%f, %f, %f)\n",s1.pos.x, s1.pos.y, s1.r);
    printf("\nmini3 (%f, %f, %f)\n",s.pos.x, s.pos.y, s.r);

    doc << svg::Circle(svg::Point(s.pos.x, s.pos.y), s.r*2, svg::Fill(), svg::Stroke(2, svg::Color(0, 0, 0)));

    for(int i = 0; i < v.size(); ++i)
        doc << svg::Circle(svg::Point(v[i].x, v[i].y), 6, svg::Fill(svg::Color(255, 0, 0)), svg::Stroke());

    doc.save();
    return 0;
}
