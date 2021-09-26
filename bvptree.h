#pragma once

#include <algorithm>
#include <random>
#include <stack>
#include <cmath>

#include "geo.h"
#include "lambert.h"

namespace geo
{
    struct BregmanKL
    {
        float operator()(const geo::Point& a, const geo::Point & b) const
        {
            return a.x*std::log(a.x/b.x)+a.y*std::log(a.y/b.y)+b.x-a.x+b.y-a.y;
        }
    };

    // cf https://stackoverflow.com/questions/4313992/methods-for-implementing-contour-plotting
    template<typename Divergence>
    std::vector<Point> ball(const geo::Point& p, float tau, Divergence d)
    {
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
                float fa = d(a,p);
                float fb = d(b,p);

                if( (fa-tau) * (fb-tau) < 0 )
                {
                    float t = fa-tau > 0 ? (fa-tau) / (fa-fb) : (fb-tau) / (fb-fa);
                    v.push_back(lerp(a,b,t));
                }
            }
            // vertical
            {
                geo::Point a = geo::Point(dx0,dy0);
                geo::Point b = geo::Point(dx0,dy1);
                float fa = d(a,p);
                float fb = d(b,p);

                if( (fa-tau) * (fb-tau) < 0 )
                {
                    float t = fa-tau > 0 ? (fa-tau) / (fa-fb) : (fb-tau) / (fb-fa);
                    v.push_back(lerp(a,b,t));
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

        auto compare = [&p0](const geo::Point& a, const geo::Point& b) 
        {
            float val = angle(p0, a, b);
            if( val < 0 ) return true;
            if( val > 0 ) return false;
            return (distance2(p0, b) >= distance2(p0, a));
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

    std::vector<geo::Point> klball(const geo::Point& p, float tau)
    {
        int nbsteps=250;
        std::vector<geo::Point> output(nbsteps*4);
        float du=tau/(float)nbsteps;

        // parametric: x1, y1, x2, y2
            // x1=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
            // y1=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
            // x2=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
            // y2=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
        
        // top left quadrant: (x1, y2)
        for (int i=0; i<=nbsteps; ++i)
        {
            float u = tau-i * du;
            float x=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
            output[i] = geo::Point(x,y);
        }
        // top right quadrant: (x2, y2)
        for (int i=0; i<=nbsteps; ++i)
        {
            float u = i * du;
            float x=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertWm1(-std::exp(-(tau-u)/p.y-1));
            output[nbsteps+i] = geo::Point(x,y);
        }
        // bottom right quadrant: (x2, y1)
        for (int i=0; i<=nbsteps; ++i)
        {
            float u = tau-i * du;
            float x=-p.x*func::LambertWm1(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
            output[2*nbsteps+i] = geo::Point(x,y);
        }
        // bottom left quadrant: (x1, y1)
        for (int i=0; i<=nbsteps; ++i)
        {
            float u = i * du;
            float x=-p.x*func::LambertW0(-std::exp(-u/p.x-1));
            float y=-p.y*func::LambertW0(-std::exp(-(tau-u)/p.y-1));
            output[3*nbsteps+i] = geo::Point(x,y);
        }
        
        return output;
    }
}

// cf: https://hal.archives-ouvertes.fr/hal-00481723/document
namespace bvptree 
{
    struct Node 
    {
        Node(int id, double r = 0., int left = -1, int right = -1)
        : id(id), r(r), left(left), right(right) {}
        int id;
        double r;
        int left;
        int right;
    };

    class Tree
    {
    public:
        
        static float bregmanKL(const geo::Point& a, const geo::Point & b)
        {
            return a.x*std::log(a.x/b.x)+a.y*std::log(a.y/b.y)+b.x-a.x+b.y-a.y;
        }

        std::vector<geo::Point> points;
        std::vector<Node> nodes;
        int root;
        std::mt19937 rng;

        Tree(){}

        Tree( const std::vector<geo::Point>& data ) : points(data)
        {
            root = build(0, points.size());
        }

        int build(int begin, int end) 
        {
            if( begin >= end ) return -1;

            if(begin + 1 == end) 
            {
                nodes.push_back(Node(begin));
                return nodes.size() - 1;
            }

            std::uniform_int_distribution<int> uni(begin, end - 1);
            std::swap(points[begin], points[uni(rng)]);

            geo::Point p = points[begin];
            auto compare = [&p](const geo::Point& a, const geo::Point& b) 
            {
                return bregmanKL(p, a) < bregmanKL(p, b);
            };

            int mid = (end + begin) / 2;
            std::nth_element( points.data() + begin+1, points.data() + mid, points.data() + end, compare);

            double r = bregmanKL(points[begin], points[mid]);          
            int left = build(begin+1, mid);
            int right = build(mid, end);  

            nodes.push_back( Node(begin, r, left, right) );
            return (int) nodes.size() - 1;
        }
    };
}