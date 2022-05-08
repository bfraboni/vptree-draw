#pragma once

#include <algorithm>
#include <random>
#include <stack>
#include <cmath>

#include "geo.h"
#include "lambert.h"

// new Bregman KL ball parametric contour
// solve for x given the Bregman ball (c, tau), B_KL(c,x) = tau 
// using the Lambert W function
namespace geo
{
    // Frank Nielsen's parametric version of KL-balls
    std::vector<geo::Point> klball(const geo::Point& p, const float tau, const int nb = 250)
    {
        std::vector<geo::Point> output(nb*4);
        float du=tau/(float)(nb-1);

        // parametric: x1, y1, x2, y2
        // x1=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
        // y1=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
        // x2=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
        // y2=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
        
        // top left quadrant: (x1, y2)
        for (int i=0; i<nb; ++i)
        {
            float u = tau-(i+0.5f) * du;
            float x=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
            output[i] = geo::Point(x,y);
        }
        // top right quadrant: (x2, y2)
        for (int i=0; i<nb; ++i)
        {
            float u = (i+0.5f) * du;
            float x=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
            output[nb+i] = geo::Point(x,y);
        }
        // bottom right quadrant: (x2, y1)
        for (int i=0; i<nb; ++i)
        {
            float u = tau-(i+0.5f) * du;
            float x=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
            output[2*nb+i] = geo::Point(x,y);
        }
        // bottom left quadrant: (x1, y1)
        for (int i=0; i<nb; ++i)
        {
            float u = (i+0.5f) * du;
            float x=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
            output[3*nb+i] = geo::Point(x,y);
        }
        
        return output;
    }

    // Frank Nielsen's parametric version of IS-balls
    std::vector<geo::Point> isball(const geo::Point& p, const float tau, const int nb = 250)
    {
        std::vector<geo::Point> output(nb*4);
        float du=tau/(float)(nb-1);

        // parametric: x1, y1, x2, y2
        // x1=-p.x/func::LambertWm1(-std::exp(-u-1));
        // y1=-p.y/func::LambertWm1(-std::exp(-(tau-u)-1));
        // x2=-p.x/func::LambertW0(-std::exp(-u-1));
        // y2=-p.y/func::LambertW0(-std::exp(-(tau-u)-1));
        
        // top left quadrant: (x1, y2)
        for (int i=0; i<nb; ++i)
        {
            float u = tau-(i+0.5) * du;
            float x=-p.x/func::LambertWm1(-std::exp(-u-1));
            float y=-p.y/func::LambertW0(-std::exp(-(tau-u)-1));
            output[i] = geo::Point(x,y);
        }
        // top right quadrant: (x2, y2)
        for (int i=0; i<nb; ++i)
        {
            float u = (i+0.5) * du;
            float x=-p.x/func::LambertW0(-std::exp(-u-1));
            float y=-p.y/func::LambertW0(-std::exp(-(tau-u)-1));
            output[nb+i] = geo::Point(x,y);
        }
        // bottom right quadrant: (x2, y1)
        for (int i=0; i<nb; ++i)
        {
            float u = tau-(i+0.5) * du;
            float x=-p.x/func::LambertW0(-std::exp(-u-1));
            float y=-p.y/func::LambertWm1(-std::exp(-(tau-u)-1));
            output[2*nb+i] = geo::Point(x,y);
        }
        // bottom left quadrant: (x1, y1)
        for (int i=0; i<nb; ++i)
        {
            float u = (i+0.5) * du;
            float x=-p.x/func::LambertWm1(-std::exp(-u-1));
            float y=-p.y/func::LambertWm1(-std::exp(-(tau-u)-1));
            output[3*nb+i] = geo::Point(x,y);
        }
        
        return output;
    }
}

// old "handmade" Bregman ball parametric contour
// 1- grid marching + bisection + interpolation to get unordered points on the contour
// 2- convex hull construction to get the final contour
// expensive computations
namespace geo
{        
    static float bregmanKL(const geo::Point& a, const geo::Point & b)
    {
        return a.x*std::log(a.x/b.x)+a.y*std::log(a.y/b.y)+b.x-a.x+b.y-a.y;
    }
    
    struct BregmanKL
    {
        float operator()(const geo::Point& a, const geo::Point & b) const
        {
            return bregmanKL(a,b);
        }
    };

    static float bregmanIS(const geo::Point& a, const geo::Point & b)
    {
        return a.x/b.x-std::log(a.x/b.x)+a.y/b.y-std::log(a.y/b.y)-2;
    }
    
    struct BregmanIS
    {
        float operator()(const geo::Point& a, const geo::Point & b) const
        {
            return bregmanIS(a,b);
        }
    };

    // cf https://stackoverflow.com/questions/4313992/methods-for-implementing-contour-plotting
    // - grid marching + bisection + interpolation 
    // - largely more robust with bisection method + final interpolation than interpolation only
    // - care with the bregman KL non symmetry
    template<typename Divergence = geo::BregmanKL>
    std::vector<Point> ball(const geo::Point& p, float tau)
    {
        Divergence d;
        int grid = 4000;
        float size = std::max(100.f,1000 * std::log(tau));
        std::vector<geo::Point> v;
        for(int i = -grid/2; i < grid/2; i++)
        for(int j = -grid/2; j < grid/2; j++)
        {
            float dx0 = p.x + float(i) * size / grid; 
            float dx1 = p.x + float(i+1) * size / grid; 
            float dy0 = p.y + float(j) * size / grid; 
            float dy1 = p.y + float(j+1) * size / grid; 

            if( dx0 < 0 || dy0 < 0 ) continue;
            if( dx1 < 0 || dy1 < 0 ) continue;
            // horizontal
            {
                geo::Point a = geo::Point(dx0,dy0);
                geo::Point b = geo::Point(dx1,dy0);
                // Warning KL divergence is not symmetric
                // should pass the center as first parameter
                float fa = d(p,a);
                float fb = d(p,b);
                float da = fa-tau;
                float db = fb-tau;

                if( da * db < 0 )
                {
                    // find the point by bisection
                    while(da > 1e-6 && db > 1e-6)
                    {
                        geo::Point x = (a+b)/2;
                        float fx = d(p,x);
                        float dx = fx-tau;
                        if( dx >= 0 && da >= 0 )
                        {
                            a = x;
                            fa = fx;
                            da = dx;
                        }
                        else
                        {
                            b = x;
                            fb = fx;
                            db = dx;
                        }
                    }
                    // v.push_back(a);

                    float t = fa-tau > 0 ? da / (fa-fb) : db / (fb-fa);
                    geo::Point x = fa-tau > 0 ? lerp(a,b,t) : lerp(b,a,t);
                    v.push_back(x);
                }
            }
            // vertical
            {
                geo::Point a = geo::Point(dx0,dy0);
                geo::Point b = geo::Point(dx0,dy1);
                // Warning KL divergence is not symmetric
                // should pass the center as first parameter
                float fa = d(p,a);
                float fb = d(p,b);
                float da = fa-tau;
                float db = fb-tau;

                if( da * db < 0 )
                {
                    // find the point by bisection
                    while(da > 1e-6 && db > 1e-6)
                    {
                        geo::Point x = (a+b)/2;
                        float fx = d(p,x);
                        float dx = fx-tau;
                        if( dx >= 0 && da >= 0 )
                        {
                            a = x;
                            fa = fx;
                            da = dx;
                        }
                        else
                        {
                            b = x;
                            fb = fx;
                            db = dx;
                        }
                    }
                    // v.push_back(a);

                    float t = fa-tau > 0 ? da / (fa-fb) : db / (fb-fa);
                    geo::Point x = fa-tau > 0 ? lerp(a,b,t) : lerp(b,a,t);
                    v.push_back(x);
                }
            }
        }

        return v;
    }

    float angle(const geo::Point& a, const geo::Point& b, const geo::Point& c)
    {
        // colinear = 0, clock > 0 or counterclock wise < 0
        return (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y); 
    }

    // cf: https://en.wikipedia.org/wiki/Introduction_to_Algorithms
    geo::Point next_to_top(std::stack<geo::Point>& stack)
    {
        geo::Point top = stack.top();
        stack.pop();
        geo::Point second = stack.top();
        stack.push(top);
        return second;
    }

    // https://en.wikipedia.org/wiki/Graham_scan
    std::vector<Point> hull(const std::vector<Point>& points)
    {
        // copy data - ugly to change later
        std::vector<geo::Point> v(points);
        if( points.empty() ) return v;
        
        int ymin = v[0].y, min = 0;
        for (int i = 1; i < v.size(); ++i)
            if((v[i].y < ymin) || (ymin == v[i].y && v[i].x < v[min].x))
                ymin = v[i].y, min = i;

        std::swap(v[0], v[min]);
        geo::Point p0 = v[0];

        BregmanKL breg;
        auto compare = [&p0, &breg](const geo::Point& a, const geo::Point& b) 
        {
            float val = angle(p0, a, b);
            if( val < 0 ) return true;
            if( val > 0 ) return false;
            return (breg(p0, b) >= breg(p0, a));
            // return (distance2(p0, b) >= distance2(p0, a));
        };

        std::sort(v.begin(), v.end(), compare);

        std::stack<geo::Point> stack({v[0],v[1],v[2]});
        for (int i = 3; i < v.size(); i++)
        {
            while(stack.size() > 1 && angle(next_to_top(stack), stack.top(), v[i]) >= 0)
                stack.pop();
            stack.push(v[i]);
        }

        v.clear();
        while(!stack.empty())
        {
            v.push_back(stack.top());
            stack.pop();
        }
        return v;
    }
}