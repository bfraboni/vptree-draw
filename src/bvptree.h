#pragma once

#include <algorithm>
#include <random>
#include <stack>
#include <cmath>

// Bregman KL balls utility functions 
#include "bregman.h"

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

    template<typename Divergence = geo::BregmanKL>
    class Tree
    {
    public:
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
            Divergence div;
            auto compare = [&p, &div](const geo::Point& a, const geo::Point& b) 
            {
                return div(p, a) < div(p, b);
            };

            int mid = (end + begin) / 2;
            std::nth_element( points.data() + begin+1, points.data() + mid, points.data() + end, compare);

            double r = div(p, points[mid]);          
            int left = build(begin+1, mid);
            int right = build(mid, end);  

            nodes.push_back( Node(begin, r, left, right) );
            return (int) nodes.size() - 1;
        }
    };
}