#pragma once

#include "geo.h"

namespace quadtree
{
    struct Node
    {
        int ne = -1, nw = -1, sw = -1, se = -1, id = -1;
        geo::Box bounds;
        Node(const geo::Box& bounds) : bounds(bounds){}
        bool leaf() const { return id < 0; }
    };

    class Tree
    {
    public:
        std::vector<geo::Point> points;
        std::vector<Node> nodes;
        int root;

        Tree(){}

        Tree( const std::vector<geo::Point>& data, const geo::Box& bounds ) : points(data)
        {
            // init root node
            root = 0;
            nodes.push_back(Node(bounds));
            // insert points
            #pragma omp parallel for
            for(unsigned int i = 0; i < points.size(); ++i)
                insert(root, i);
        }

        void subdivide(int node)
        {
            geo::Point pmin = nodes[node].bounds.pmin;
            geo::Point pmax = nodes[node].bounds.pmax;
            geo::Point pmid = (pmin+pmax)*0.5f;
            // north east
            nodes.push_back(Node(geo::Box(pmid, pmax)));
            nodes[node].ne = nodes.size()-1;
            // north west
            nodes.push_back(Node(geo::Box(geo::Point(pmin.x, pmid.y), geo::Point(pmid.x, pmax.y))));
            nodes[node].nw = nodes.size()-1;
            // south west
            nodes.push_back(Node(geo::Box(pmin, pmid)));
            nodes[node].sw = nodes.size()-1;
            // south east
            nodes.push_back(Node(geo::Box(geo::Point(pmid.x, pmin.y), geo::Point(pmax.x, pmid.y))));
            nodes[node].se = nodes.size()-1;
        }

        bool insert(int node, int id)
        {
            const geo::Point& p = points[id]; 

            if( !nodes[node].bounds.contains(p) ) return false;

            if( nodes[node].id < 0 )
            {
                nodes[node].id = id;
                return true;
            }

            if(nodes[node].nw < 0) subdivide(node);
            if(insert(nodes[node].ne, id)) return true;
            if(insert(nodes[node].nw, id)) return true;
            if(insert(nodes[node].sw, id)) return true;
            if(insert(nodes[node].se, id)) return true;
            return false;   
        }
    };
}