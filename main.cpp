#include <iostream> 
#include <random>
#include <cassert>

#include "vptree.hpp" 
#include "simple_svg_1.0.0.hpp" 
#include "cavc/polylinecombine.hpp"

using namespace svg;

namespace svg
{
    class Arc : public Shape
    {
    public:
        Arc(Point const & center, double diameter, double start, double end, 
            Fill const & fill, Stroke const & stroke = Stroke())
            : Shape(fill, stroke), center(center), radius(diameter / 2), start(start), end(end) { }

        // arc example from : https://www.w3.org/TR/SVG/paths.html#PathDataEllipticalArcCommands 
        // <path d="M 125,75 a100,50 0 ?,? 100,50" style="fill:none; stroke:red; stroke-width:6"/>
        std::string toString(Layout const & layout) const
        {
            double sr = (start-90.0) * M_PI / 180.0;
            Point ps(center.x + radius * std::cos(sr), center.y + radius * std::sin(sr));
            
            double se = (end-90.0) * M_PI / 180.0;
            Point pe(center.x + radius * std::cos(se), center.y + radius * std::sin(se));
            
            std::stringstream ss;
            ss << elemStart("path");
            ss << "d=\"";
            ss << "M";
            ss << translateX(ps.x, layout) << "," << translateY(ps.y, layout);
            ss << "A";
            ss << translateScale(radius, layout) << "," << translateScale(radius, layout);
            ss << " 0 ";
            ss << (end - start <= 180.0 ? "0," : "1,");
            ss << "0 ";
            ss << translateX(pe.x, layout) << "," << translateY(pe.y, layout);
            ss << "\" ";
            ss << fill.toString(layout) << stroke.toString(layout) << emptyElemEnd();
            return ss.str();
        }
        void offset(Point const & offset)
        {
            center.x += offset.x;
            center.y += offset.y;
        }
    private:
        Point center;
        double radius, start, end;
    };

    class CavcPoly : public Shape
    {
    public:
        CavcPoly(const cavc::Polyline<double>& poly, Fill const & fill, Stroke const & stroke = Stroke())
            : Shape(fill, stroke), poly(poly) { }

        std::string toString(Layout const & layout) const
        {
            std::stringstream ss;
            for( int i = 0; i < poly.size(); i++)
            {
                const auto& p1 = poly.vertexes()[i];
                const auto& p2 = poly.vertexes()[(i+1)%poly.size()];

                if(!p1.bulgeIsZero())
                {
                    auto arc = arcRadiusAndCenter(p1, p2);

                    double start = (angle(arc.center, p1.pos()) + M_PI) * 180 / M_PI;
                    double end = (angle(arc.center, p2.pos()) + M_PI) * 180 / M_PI;
                    double radius = arc.radius;
                    Point ps(p1.pos().x(), p1.pos().y());
                    Point pe(p2.pos().x(), p2.pos().y());

                    ss << elemStart("path");
                    ss << "d=\"";
                    ss << "M";
                    ss << translateX(ps.x, layout) << "," << translateY(ps.y, layout);
                    ss << "A";
                    ss << translateScale(radius, layout) << "," << translateScale(radius, layout);
                    ss << " 0 ";
                    // ss << (end - start <= 180.0 ? "1," : "0,");
                    ss << "0,";
                    // ss << "1 ";
                    ss << (p1.bulge() > 0 ? "0," : "1,");;
                    ss << translateX(pe.x, layout) << "," << translateY(pe.y, layout);
                    ss << "\" ";
                    ss << fill.toString(layout) << stroke.toString(layout) << emptyElemEnd();
                }
            }

            return ss.str();
        }
        void offset(Point const & offset) { }
    private:
        const cavc::Polyline<double>& poly;
    };
}

void draw(const vpt::VpTree& vptree, svg::Document& doc)
{
    for(int i = 0; i < vptree.items_.size(); i++)
    {
        const auto& point = vptree.items_[i].first;
        doc << Circle(Point(point[0], point[1]), 1, Fill(Color(0, 0, 0)), Stroke());

        // circle border
        int id = vptree.items_[i].second;
        const auto& node = vptree.nodes_[id];
        doc << Circle(Point(point[0], point[1]), 2*node.threshold, Fill(), Stroke(1, Color(125, 125, 125)));

        // doc << Circle(Point(point[0], point[1]), 1, Fill(Color(100, 200, 120)), Stroke(1, Color(200, 250, 150)));
    }
}   

/*
void drawVPTree(int node, const cavc::Polyline<double>& path, const svg::Color& color) 
{
    if(node == vpt::Node::Leaf)
        return;

    auto inside = area;
    auto outside = area;

    if( node == 0 ) 
    {
        // construct new circle
        const auto& n = vptree.nodes_[node];
        double radius = n.threshold;
        double cx =  vptree.items_[n.item][0];
        double cy =  vptree.items_[n.item][1];

        // closed polyline representing a circle
        cavc::Polyline<double> circle;
        circle.addVertex(cx-radius, cy, radius);
        circle.addVertex(cx+radius, cy, radius);
        circle.isClosed() = true;

        // draw circle

        // compute intersection with enclosing path
        node.circle = circle;
        circle.strokeColor = color;
        circle.strokeWidth = strokeWidth;

        trash.push(circle);

        inside = path.intersect(circle);
        outside = path.subtract(circle);


//     cavc::CombineResult<double> unionResult = combinePolylines(circle, rectangle, PlineCombineMode::Union);
//     cavc::CombineResult<double> excludeResult = combinePolylines(circle, rectangle, PlineCombineMode::Exclude);
//     cavc::CombineResult<double> intersectResult = combinePolylines(circle, rectangle, PlineCombineMode::Intersect);
    } 
    else 
    {
        var circle = new vppaper.Path.Circle(node.vp, node.mu);
        var intersect = area.intersect(circle);
        intersect.strokeColor = color;
        intersect.strokeWidth = strokeWidth;

        trash.push(intersect);
        trash.push(circle);

        inside = intersect;
        outside = area.subtract(intersect);
    }

    color.hue += 10;
    if (drawType != "all") 
    {
        vpsteps.push([node.left, inside, color]);
        vpsteps.push([node.right, outside, color]);
    } 
    else 
    {
        drawVPTree(node.left, inside, color);
        drawVPTree(node.right, outside, color);    
    }
}
*/

// using namespace cavc;

// int main(int argc, char *argv[]) 
// {
//     (void)argc;
//     (void)argv;

//     // closed polyline representing a circle
//     cavc::Polyline<double> circle;
//     circle.addVertex(0, 1, 1);
//     circle.addVertex(10, 1, 1);
//     circle.isClosed() = true;

//     // closed polyline representing a rectangle (overlaps with the circle)
//     cavc::Polyline<double> rectangle;
//     rectangle.addVertex(3, -10, 0);
//     rectangle.addVertex(6, -10, 0);
//     rectangle.addVertex(6, 10, 0);
//     rectangle.addVertex(3, 10, 0);
//     rectangle.isClosed() = true;

//     cavc::CombineResult<double> unionResult = combinePolylines(circle, rectangle, PlineCombineMode::Union);
//     cavc::CombineResult<double> excludeResult = combinePolylines(circle, rectangle, PlineCombineMode::Exclude);
//     cavc::CombineResult<double> intersectResult = combinePolylines(circle, rectangle, PlineCombineMode::Intersect);
//     cavc::CombineResult<double> xorResult = combinePolylines(circle, rectangle, PlineCombineMode::XOR);

//     return 0;
// }

int main()
{   
    srand(1234);

    int dx = 1000, dy = 1000;
    int scale = 1;
    Dimensions dimensions(dx, dy);
    Layout layout(dimensions, Layout::BottomLeft, scale, Point(0, 0));
    Document doc("vptree.svg", layout);

    // background.
    Polygon border(Fill(Color(255, 255, 255)), Stroke());
    border << Point(0, 0) << Point(dimensions.width, 0) << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    doc << border;

    // auto arc = Arc(Point(dx/2, dy/2), 250, 0, 90, Fill(Color(255, 0, 0)), Stroke(5, Color(0, 255, 0))) ;
    // auto arc = Arc(Point(dx/2, dy/2), 250, 0, 270, Fill(), Stroke(5, Color(0, 255, 0))) ;
    // std::cout << arc.toString(layout) << std::endl;
    // doc << arc;

    // init random points
    int nb = 1000;
    auto points = std::vector<std::vector<double>>(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = drand48() * dx/2 + dx/4;
        double y = drand48() * dy/2 + dy/4;
        points[i] = std::vector<double>({x,y,0});
    }

    // create a vantage points tree
    vpt::VpTree t1(points);

    // draw vptree to svg
    // draw(t1, doc);

    // test draw path intersection

    //doc << Circle(Point(dx/2, dy/2), dx, Fill(), Stroke(10, Color(255, 0, 255)));
    auto n = t1.root();
    auto p = t1.items_[n.item].first;
    double diam = 2 * n.threshold;
    //doc << Circle(Point(p[0], p[1]), diam, Fill(), Stroke(10, Color(0, 255, 0))); 
    //doc << Circle(Point(p[0], p[1]), 10, Fill(Color(0, 255, 0)), Stroke());
    

    // closed polyline representing a circle
    cavc::Polyline<double> circle1;
    circle1.addVertex(150-150, 150, 1);
    circle1.addVertex(150+150, 150, 1);
    circle1.isClosed() = true;

    cavc::Polyline<double> circle2;
    circle2.addVertex(0, dy/2, 1);
    circle2.addVertex(dx, dy/2, 1);
    circle2.isClosed() = true;

    cavc::CombineResult<double> excludeResult = combinePolylines(circle2, circle1, cavc::PlineCombineMode::Exclude);
    cavc::CombineResult<double> intersectResult = combinePolylines(circle2, circle1, cavc::PlineCombineMode::Intersect);

    std::cout << "intersect: " << intersectResult.remaining.size() << std::endl;
    std::cout << " remains " << intersectResult.subtracted.size() << " holes " << std::endl;
    std::cout << "exclude: " << excludeResult.remaining.size() << std::endl;
    std::cout << " remains " << excludeResult.subtracted.size() << " holes " << std::endl;

    for(const auto& pl : intersectResult.remaining)
    {
        doc << svg::CavcPoly(pl, Fill(), Stroke(2, Color(0, 0, 0)));
    }
    // for(const auto& pl : intersectResult.subtracted)
    // {
    //     doc << svg::CavcPoly(pl, Fill( Color(100, 100, 100) ), Stroke(2, Color(0, 0, 0)));
    // }

    for(const auto& pl : excludeResult.remaining)
    {
        auto svgpl = svg::CavcPoly(pl, Fill(), Stroke(2, Color(0, 0, 0)));
        doc << svgpl;
        std::cout << svgpl.toString(layout) << std::endl;
    }

    for(const auto& pl : excludeResult.subtracted)
    {
        auto svgpl = svg::CavcPoly(pl, Fill(), Stroke(2, Color(0, 0, 255)));
        doc << svgpl;
        std::cout << svgpl.toString(layout) << std::endl;
    }

    // n = t1.nodes_[t1.nodes_[0].left];
    // radius = 2 * n.threshold;
    // cx =  t1.items_[n.item].first[0];
    // cy =  t1.items_[n.item].first[1];
    // circle = Circle(Point(cx, cy), radius, Fill(), Stroke(5, Color(255, 0, 0)));
    // doc << circle;
    // doc << Circle(Point(cx, cy), 4, Fill(Color(255, 0, 0)), Stroke());
    
    // n = t1.nodes_[t1.nodes_[0].right];
    // radius = 2 * n.threshold;
    // cx =  t1.items_[n.item].first[0];
    // cy =  t1.items_[n.item].first[1];
    // circle = Circle(Point(cx, cy), radius, Fill(), Stroke(5, Color(0, 0, 255)));
    // doc << circle;
    // doc << Circle(Point(cx, cy), 4, Fill(Color(0, 0, 255)), Stroke());

    // save svg
    doc.save();

    return 0;
}
