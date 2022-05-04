#pragma once

#include <algorithm>

#include "geo.h"

namespace kdtree
{
    struct Node
    {   
        int id, left, right;
        Node( int id, int left, int right ) : id(id), left(left), right(right) {}
    };

    struct Tree
    {
        std::vector<Node> nodes;
        std::vector<geo::Point> points;
        int root = -1;

        Tree( const std::vector<geo::Point>& points ) : points(points)
        {
            root = build(0, points.size(), 0);
        }

        int build( const int begin, const int end, const int depth )
        {
            if( begin >= end ) return -1;

            int axis = depth % 2;
            int len = end - begin;
            int mid = begin + len / 2;

            auto compare = [axis](const geo::Point& a, const geo::Point& b) 
            {
                return a[axis] < b[axis];
            };

            std::nth_element(points.data() + begin, points.data() + mid, points.data() + end, compare);

            int left = build( begin, mid, depth + 1 );
            int right = build( mid + 1, end, depth + 1 );

            nodes.push_back( Node( mid, left, right ) );
            
            return (int) nodes.size() - 1;
        }
    };
}