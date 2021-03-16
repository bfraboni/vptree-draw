#pragma once

#include "bvhbox.h"
#include "bvhsphere.h"
#include "kdtree.h"
#include "quadtree.h"
#include "vptree.h"

// svg drawing library
#include "simple_svg_1.0.0.hpp" 

// polyline contour intersections
#include "cavc/polylinecombine.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
// polyline contour svg drawing
#include "simple_svg_extend.h" 

// utility 
svg::Color rotate( const svg::Color &in, double hue)
{
    double U = cos(std::fmod(hue,360.)*M_PI/180.);
    double W = sin(std::fmod(hue,360.)*M_PI/180.);
    const double rgb[3] = {in.red/255., in.green/255., in.blue/255.};
    double res[3] = {
        ((.299+.701*U+.168*W)*rgb[0] + (.587-.587*U+.330*W)*rgb[1] + (.114-.114*U-.497*W)*rgb[2])*255.,
        ((.299-.299*U-.328*W)*rgb[0] + (.587+.413*U+.035*W)*rgb[1] + (.114-.114*U+.292*W)*rgb[2])*255.,
        ((.299-.3*U+1.25*W)*rgb[0] + (.587-.588*U-1.05*W)*rgb[1] + (.114+.886*U-.203*W)*rgb[2])*255.
    };
    return svg::Color(res[0], res[1], res[2]);
}

namespace bvhsphere
{
    void draw( const bvhsphere::Tree& tree, int depth, int node, svg::Document& doc)
    {
        const auto& n = tree.nodes[node];
        
        // hue rotated color 
        svg::Color color = rotate(svg::Color(255,120,80), depth*30);

        // if leaf
        if(n.leaf())
        {
            const geo::Point& p = tree.points[n.leaf_begin()]; 
            doc << svg::Circle(svg::Point(p.x, p.y), 4, svg::Fill(color), svg::Stroke());

            // std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
            // std::cout << "circle: " << p.x << " " << p.y << " " << 4 << std::endl;
        }
        // if node
        else
        {
            const geo::Point& p = n.bounds.pos;
            const float radius = n.bounds.r; 
            const float diameter = 2 * radius; 

            // std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
            // std::cout << "circle: " << p.x << " " << p.y << " " << radius << std::endl;

            // draw circle center
            doc << svg::Circle(svg::Point(p.x, p.y), diameter, svg::Fill(), svg::Stroke(2, color));

            // draw left subtree
            draw(tree, depth+1, n.left, doc);

            // draw right subtree
            draw(tree, depth+1, n.right, doc);
        }
    }
}

namespace bvhbox
{
    void draw(const bvhbox::Tree& tree, int depth, int node, svg::Document& doc)
    {
        const auto& n = tree.nodes[node];
        
        // hue rotated color 
        svg::Color color = rotate(svg::Color(255,120,80), depth*30);

        // if leaf
        if(n.leaf())
        {
            const geo::Point& p = tree.points[n.leaf_begin()]; 
            doc << svg::Circle(svg::Point(p.x, p.y), 4, svg::Fill(color), svg::Stroke());

            // std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
            // std::cout << "circle: " << p.x << " " << p.y << " " << 4 << std::endl;
        }
        // if node
        else
        {
            const geo::Point& pmin = n.bounds.pmin;
            const geo::Point& pmax = n.bounds.pmax;

            // std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
            // std::cout << "box: " << pmin << " " << pmax << " " << pmax - pmin << std::endl;

            // draw box
            svg::Polygon poly(svg::Fill(), svg::Stroke(2, color));
            poly << svg::Point(pmin.x, pmin.y) << svg::Point(pmin.x, pmax.y) << svg::Point(pmax.x, pmax.y) << svg::Point(pmax.x, pmin.y); 
            doc << poly;

            // draw left subtree
            draw(tree, depth+1, n.left, doc);

            // draw right subtree
            draw(tree, depth+1, n.right, doc);
        }
    }
}

namespace kdtree
{
    void draw(const kdtree::Tree& tree, int depth, int node, geo::Point pmin, geo::Point pmax, svg::Document& doc)
    {
        if( node < 0 ) return;

        const auto& n = tree.nodes[node];
        
        // hue rotated color 
        svg::Color color = rotate(svg::Color(255,120,80), depth*30);

        const auto& p = tree.points[n.id];
        doc << svg::Circle(svg::Point(p.x, p.y), 4, svg::Fill(color), svg::Stroke());

        int axis = depth % 2;
        int oaxis = (axis+1)%2;
        geo::Point p0, p1;
        float value = std::min(pmax[axis], std::max(pmin[axis], p[axis]));
        p0[axis] = value;
        p0[oaxis] = pmin[oaxis];
        p1[axis] = value;
        p1[oaxis] = pmax[oaxis];

        doc << svg::Line(svg::Point(p0.x, p0.y), svg::Point(p1.x, p1.y), svg::Stroke(2, color));

        draw(tree, depth+1, n.left, pmin, p1, doc);
        draw(tree, depth+1, n.right, p0, pmax, doc);
    }
}

namespace quadtree
{
    void draw(const quadtree::Tree& tree, int depth, int node, svg::Document& doc)
    {
        if( node < 0 ) return;

        const auto& n = tree.nodes[node];
        
        // hue rotated color 
        svg::Color color = rotate(svg::Color(255,120,80), depth*30);

        // std::cout << "depth: " << depth << " " << node << " " << n.ne << " " << n.nw << " " << n.sw << " " << n.se << std::endl;
        if( !n.leaf() )
        {
            const auto& p = tree.points[n.id];
            doc << svg::Circle(svg::Point(p.x, p.y), 4, svg::Fill(color), svg::Stroke());

            // std::cout << "circle: " << p.x << " " << p.y << std::endl;

            // draw north east subtree
            draw(tree, depth+1, n.ne, doc);
            // draw north west subtree
            draw(tree, depth+1, n.nw, doc);
            // draw south west subtree
            draw(tree, depth+1, n.sw, doc);
            // draw south east subtree
            draw(tree, depth+1, n.se, doc);
        }

        // draw box
        const auto& pmin = n.bounds.pmin;
        const auto& pmax = n.bounds.pmax;
        // std::cout << "poly: " << pmin << " " << pmax << std::endl;
        svg::Polygon poly(svg::Fill(), svg::Stroke(2, color));
        poly << svg::Point(pmin.x, pmin.y) << svg::Point(pmin.x, pmax.y) << svg::Point(pmax.x, pmax.y) << svg::Point(pmax.x, pmin.y); 
        doc << poly;
    }
}

namespace vptree
{
    struct Cell
    {
        std::vector<cavc::Polyline<double>> shape;
        std::vector<cavc::Polyline<double>> holes;
    };

    void draw(
        const vptree::Tree& tree, 
        int depth, 
        int node, 
        const Cell& cell,
        std::vector<svg::CavcPoly::Edge> &edgeBuffer, 
        svg::Document& doc
    )
    {
        if(node < 0) return;

        // get cell vantage circle info
        const auto& n = tree.nodes[node];
        double radius = n.r;
        double cx = tree.points[n.id][0];
        double cy = tree.points[n.id][1];
        
        // hue rotated color 
        svg::Color color = rotate(svg::Color(255,120,80), depth*30);
        // std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
        // std::cout << "circle: " << cx << " " << cy << " " << radius << std::endl;

        // on a leaf
        if( radius <= 0 )
        {
            // draw vantage point
            doc << svg::Circle(svg::Point(cx, cy), 4, svg::Fill(color), svg::Stroke());
        }
        // on a node
        else
        {
            // closed polyline of the vantage circle
            cavc::Polyline<double> circle;
            circle.addVertex(cx-radius, cy, 1);
            circle.addVertex(cx+radius, cy, 1);
            circle.isClosed() = true;

            Cell inside, outside;
            // hole handling
            if( cell.holes.size() > 0 )
            {
                for(int i = 0; i < (int)cell.holes.size(); ++i)
                {
                    cavc::CombineResult<double> interHoleCircle = combinePolylines(cell.holes[i], circle, cavc::PlineCombineMode::Intersect);
                    // if hole and circle intersects the hole will be removed next
                    if( interHoleCircle.remaining.size() > 0 )
                    {
                        // intersect with enclosing path
                        for(int j = 0; j < (int)cell.shape.size(); ++j)
                        {
                            // intersect with enclosing path
                            cavc::CombineResult<double> interShapeCircle = combinePolylines(cell.shape[j], circle, cavc::PlineCombineMode::Intersect);
                            // std::cout << "hole but intersect: " << interShapeCircle.remaining.size() << " remains " << interShapeCircle.subtracted.size() << " holes " << std::endl;
                            for(auto& path : interShapeCircle.remaining)
                            {
                                // remove hole from new outlines
                                cavc::CombineResult<double> excluNewShapeHole = combinePolylines(path, cell.holes[i], cavc::PlineCombineMode::Exclude);
                                inside.shape.insert(inside.shape.end(), excluNewShapeHole.remaining.begin(), excluNewShapeHole.remaining.end());
                                inside.holes.insert(inside.holes.end(), excluNewShapeHole.subtracted.begin(), excluNewShapeHole.subtracted.end());
                            }
                        }

                        // compute exclusion with enclosing path
                        for(int j = 0; j < (int)cell.shape.size(); ++j)
                        {
                            cavc::CombineResult<double> excluShapeCircle = combinePolylines(cell.shape[j], circle, cavc::PlineCombineMode::Exclude);
                            // std::cout << "hole but exclude: " << excluShapeCircle.remaining.size() << " remains " << excluShapeCircle.subtracted.size() << " holes " << std::endl;
                            for(auto& path : excluShapeCircle.remaining)
                            {
                                // remove hole from new outlines
                                cavc::CombineResult<double> excluNewShapeHole = combinePolylines(path, cell.holes[i], cavc::PlineCombineMode::Exclude);
                                outside.shape.insert(outside.shape.end(), excluNewShapeHole.remaining.begin(), excluNewShapeHole.remaining.end());
                                outside.holes.insert(outside.holes.end(), excluNewShapeHole.subtracted.begin(), excluNewShapeHole.subtracted.end());
                            }
                        }
                    }
                    // else the hole wil remain at next level
                    else
                    {
                        outside.holes.push_back(cell.holes[i]);

                        // compute intersection with enclosing path
                        for(int j = 0; j < (int)cell.shape.size(); ++j)
                        {
                            cavc::CombineResult<double> inter = combinePolylines(cell.shape[j], circle, cavc::PlineCombineMode::Intersect);
                            // std::cout << "hole but do not intersect: " << inter.remaining.size() << " remains " << inter.subtracted.size() << " holes " << std::endl;
                            inside.shape.insert(inside.shape.end(), inter.remaining.begin(), inter.remaining.end());
                            inside.holes.insert(inside.holes.end(), inter.subtracted.begin(), inter.subtracted.end());
                        }

                        // compute exclusion with enclosing path
                        for(int j = 0; j < (int)cell.shape.size(); ++j)
                        {
                            cavc::CombineResult<double> exclu = combinePolylines(cell.shape[j], circle, cavc::PlineCombineMode::Exclude);
                            // std::cout << "hole but do not exclude: " << exclu.remaining.size() << " remains " << exclu.subtracted.size() << " holes " << std::endl;
                            outside.shape.insert(outside.shape.end(), exclu.remaining.begin(), exclu.remaining.end());
                            outside.holes.insert(outside.holes.end(), exclu.subtracted.begin(), exclu.subtracted.end());
                        }
                    }
                }
            }
            else
            {
                // compute intersection with enclosing path
                for(int i = 0; i < (int)cell.shape.size(); ++i)
                {
                    cavc::CombineResult<double> inter = combinePolylines(cell.shape[i], circle, cavc::PlineCombineMode::Intersect);
                    // std::cout << "intersect: " << inter.remaining.size() << " remains " << inter.subtracted.size() << " holes " << std::endl;
                    inside.shape.insert(inside.shape.end(), inter.remaining.begin(), inter.remaining.end());
                    inside.holes.insert(inside.holes.end(), inter.subtracted.begin(), inter.subtracted.end());
                }

                // compute exclusion with enclosing path
                for(int i = 0; i < (int)cell.shape.size(); ++i)
                {
                    cavc::CombineResult<double> exclu = combinePolylines(cell.shape[i], circle, cavc::PlineCombineMode::Exclude);
                    // std::cout << "exclude: " << exclu.remaining.size() << " remains " << exclu.subtracted.size() << " holes " << std::endl;
                    outside.shape.insert(outside.shape.end(), exclu.remaining.begin(), exclu.remaining.end());
                    outside.holes.insert(outside.holes.end(), exclu.subtracted.begin(), exclu.subtracted.end());
                }
            }

            // draw circle center
            doc << svg::Circle(svg::Point(cx, cy), 4, svg::Fill(color), svg::Stroke());

            // draw left subtree
            draw(tree, depth+1, n.left, inside, edgeBuffer, doc);

            // draw right subtree
            draw(tree, depth+1, n.right, outside, edgeBuffer, doc);

            // draw left cell (inside)
            for(int i = 0; i < (int)inside.shape.size(); ++i)
            {
                doc << svg::CavcPoly(inside.shape[i], svg::Fill(), edgeBuffer, svg::Stroke(2, color));
            }

            // draw right cell (outside)
            for(int i = 0; i < (int)outside.shape.size(); ++i)
            {
                doc << svg::CavcPoly(outside.shape[i], svg::Fill(), edgeBuffer, svg::Stroke(2, color));
            }    
        }
    }
}
