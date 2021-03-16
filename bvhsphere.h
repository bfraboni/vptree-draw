#pragma once

#include <iostream>
#include <cassert>
#include <algorithm>

#include "geo.h"

namespace bvhsphere
{
    struct Node 
    {
        geo::Sphere bounds;
        int left;
        int right;
        
        Node& make_leaf( const geo::Sphere& _bound, const int _begin, const int _end )
        {
            assert(_begin < _end);
            bounds = _bound;
            left = _begin;
            right = - (_end - _begin);
            assert(right < 0);
            return *this;
        }
        
        Node& make_node( const geo::Sphere& _bound, const int _left, const int _right )
        {
            bounds = _bound;
            left = _left;
            right = _right;
            return *this;
        }
        
        bool internal()  const { return (right >= 0); }
        bool leaf()      const { return (right < 0); }
        int leaf_begin() const { assert(leaf()); return left; }
        int leaf_n()     const { assert(leaf()); return -right; }
        int leaf_end()   const { assert(leaf()); return left - right; }
    };

    class Tree
    {
    public:
        std::vector<geo::Point> points;
        std::vector<Node> nodes;
        int root;

        Tree(){}
        
        Tree( const std::vector<geo::Point>& data ):points(data)
        {
            // build tree
            root = build(0, points.size());
        }

        int build( const int begin, const int end )
        {
            if(begin+1>=end)
            {
                Node node; 
                node.make_leaf( geo::Sphere(points[begin]), begin, end);
                nodes.push_back(node);
                // printf("leaf : id %lu c %f %f r %f p %f %f id %d\n", nodes.size()-1, node.bounds.pos.x, node.bounds.pos.y, node.bounds.r, points[begin].x, points[begin].y, begin);
                return (int) nodes.size() - 1;
            }
            
            // englobant des centres des objets
            geo::Box box;
            for(int i = begin; i < end; ++i)
                box.insert(points[i]);
            
            // axe et coupe
            geo::Vector d = box.pmax - box.pmin;
            int axis = d.x > d.y ? 0 : 1;
            float cut = box.pmin[axis] + d[axis] / 2;
            
            // predicate
            auto compare = [&axis, &cut](const geo::Point& p) 
            {
                return p[axis] <= cut;
            };
            
            // partition
            geo::Point * pmid = std::partition(points.data() + begin, points.data() + end, compare);
            int mid = std::distance(points.data(), pmid);
            
            // cas degenere, reparti arbitrairement les objets en 2
            if(mid == begin || mid == end)
                mid = (begin + end) / 2;
            assert(mid != begin);
            assert(mid != end);
            
            // construction recursive des fils
            int left = build(begin, mid);
            int right = build(mid, end);

            // construction du noeud
            Node node;
            node.make_node( geo::Sphere(nodes[left].bounds, nodes[right].bounds), left, right );
            nodes.push_back(node);
            // printf("node : id %lu c %f %f r %f\n", nodes.size()-1, node.bounds.pos.x, node.bounds.pos.y, node.bounds.r);
            return (int) nodes.size() -1;
        }
    };
}