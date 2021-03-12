#include <iostream> 
#include <random>
#include <cassert>

#include "vptree.hpp" 
#include "simple_svg_1.0.0.hpp" 
#include "cavc/polylinecombine.hpp"

using namespace svg;

// extend simple svg with arcs and cavc::Polyline<double> support
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
                    // arc if bulge != 0
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
                    // always small arcs -> large arc flag to 0
                    ss << "0,";                           
                    // sweep is 0 if bulge is positive, 1 if negative
                    ss << (p1.bulge() > 0 ? "0," : "1,"); 
                    ss << translateX(pe.x, layout) << "," << translateY(pe.y, layout);
                    ss << "\" ";
                    ss << fill.toString(layout) << stroke.toString(layout) << emptyElemEnd();
                }
                else
                {
                    // line if bulge == 0
                    Path p(fill, stroke);
                    p << svg::Point(p1.pos().x(), p1.pos().y()) << svg::Point(p2.pos().x(), p2.pos().y());
                    ss << p.toString(layout);
                }
            }

            return ss.str();
        }
        void offset(Point const & offset) { }
    private:
        const cavc::Polyline<double>& poly;
    };
}

svg::Color rotate(const svg::Color& color, int hue) 
{
    const float cosA = std::cos((hue%360) * M_PI / 180.f);
    const float sinA = std::sin((hue%360) * M_PI / 180.f);
    const float neo[3] = {
        cosA + (1 - cosA) / 3,
        (1 - cosA) / 3 - std::sqrt(1.f / 3.f) * sinA,
        (1 - cosA) / 3 + std::sqrt(1.f / 3.f) * sinA,
    };
    const int rgb[3] = {color.red, color.green, color.blue};
    const int res[3] = {
        std::min(255, std::max(0, (int)(rgb[0] * neo[0] + rgb[1] * neo[1] + rgb[2] * neo[2]))),
        std::min(255, std::max(0, (int)(rgb[0] * neo[2] + rgb[1] * neo[0] + rgb[2] * neo[1]))),
        std::min(255, std::max(0, (int)(rgb[0] * neo[1] + rgb[1] * neo[2] + rgb[2] * neo[0]))),
    };
    return svg::Color(res[0], res[1], res[2]);
}

struct Cell
{
    std::vector<cavc::Polyline<double>> shape;
    std::vector<cavc::Polyline<double>> holes;
};

void drawVPTree(const vpt::VpTree& vptree, int depth, int node, const Cell& cell, svg::Document& doc) 
{
    if(node == vpt::VpTree::Node::Leaf)
        return;

    // get cell vantage circle info
    const auto& n = vptree.nodes_[node];
    double radius = n.threshold;
    double cx = vptree.items_[n.item].first[0];
    double cy = vptree.items_[n.item].first[1];

    std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
    std::cout << "circle: " << cx << " " << cy << " " << radius << std::endl;

    // on a leaf
    if( radius <= 0 )
    {
        // draw vantage point
        auto svgpl = svg::Circle(svg::Point(cx, cy), 4, Fill(svg::Color::Defaults::Black), Stroke());
        doc << svgpl;
    }
    // on a node
    else
    {
        // closed polyline of the vantage circle
        cavc::Polyline<double> circle;
        circle.addVertex(cx-radius, cy, 1);
        circle.addVertex(cx+radius, cy, 1);
        circle.isClosed() = true;

        // compute intersection with enclosing path
        Cell inside;
        for(int i = 0; i < cell.shape.size(); ++i)
        {
            cavc::CombineResult<double> inter = combinePolylines(circle, cell.shape[i], cavc::PlineCombineMode::Intersect);
            std::cout << "intersect: " << inter.remaining.size() << " remains " << inter.subtracted.size() << " holes " << std::endl;
            if( inter.remaining.size() > 0 )
            {
                inside.shape.insert(inside.shape.end(), inter.remaining.begin(), inter.remaining.end());
            }
        }

        // compute exclusion with enclosing path
        Cell outside;
        for(int i = 0; i < cell.shape.size(); ++i)
        {
            cavc::CombineResult<double> exclu = combinePolylines(cell.shape[i], circle, cavc::PlineCombineMode::Exclude);
            std::cout << "exclude: " << exclu.remaining.size() << " remains " << exclu.subtracted.size() << " holes " << std::endl;
            if( exclu.remaining.size() > 0 )
            {
                outside.shape.insert(outside.shape.end(), exclu.remaining.begin(), exclu.remaining.end());
            }

            // if there is a hole -> TODO
            if( exclu.subtracted.size() > 0 )
            {
                std::cout << "\tFire in the hole\n" << std::endl;
                outside.holes.insert(outside.holes.end(), exclu.subtracted.begin(), exclu.subtracted.end());
            }
        }

        // draw inside cell path
        svg::Color color = rotate(svg::Color::Defaults::Red, depth * 30);
        for(int i = 0; i < inside.shape.size(); ++i)
        {
            auto svgpl = svg::CavcPoly(inside.shape[i], Fill(), Stroke(2, color));
            doc << svgpl;
        }

        // draw outside cell path
        for(int i = 0; i < outside.shape.size(); ++i)
        {
            auto svgpl = svg::CavcPoly(outside.shape[i], Fill(), Stroke(2, color));
            doc << svgpl;
        }   
        
        // draw subtrees
        drawVPTree(vptree, depth+1, n.left, inside, doc);
        drawVPTree(vptree, depth+1, n.right, outside, doc); 
    }
}

static float rand1D() 
{   
    using G = std::default_random_engine;    
    using D = std::uniform_real_distribution<float>; 
    static const std::size_t seed = 123456789;   
    static G gen(seed);
    static D dist(0.f, 1.f);
    return dist(gen);
}

int main()
{   
    int dx = 1000, dy = 1000;
    int scale = 1;
    Dimensions dimensions(dx, dy);
    Layout layout(dimensions, Layout::BottomLeft, scale, Point(0, 0));
    Document doc("vptree.svg", layout);

    // background white rectangle
    Polygon border(Fill(Color(255, 255, 255)), Stroke());
    border << Point(0, 0) << Point(dimensions.width, 0) << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    doc << border;

    // init random points
    int nb = 1000;
    auto points = std::vector<std::vector<double>>(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = rand1D() * dx;
        double y = rand1D() * dy;
        points[i] = std::vector<double>({x,y,0});
    }

    // create a vantage points tree
    vpt::VpTree vptree(points);

    // init cell representing the canvas
    cavc::Polyline<double> shape1;
    shape1.addVertex(0, 0, 0); 
    shape1.addVertex(dx-0, 0, 0); 
    shape1.addVertex(dx-0, dy-0, 0); 
    shape1.addVertex(0, dy-0, 0);
    shape1.isClosed() = true;
    
    Cell cell;
    cell.shape.push_back(shape1);
    
    // draw VP tree
    drawVPTree(vptree, 0, 0, cell, doc);

    // save svg
    doc.save();

    return 0;
}
