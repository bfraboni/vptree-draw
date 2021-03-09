#include <iostream> 
#include <random>

#include "vptree.hpp" 
#include "simple_svg_1.0.0.hpp" 

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
}

void draw(const vpt::VpTree& vptree, svg::Document& doc)
{
    for(int i = 0; i < vptree.items_.size(); i++)
    {
        const auto& point = vptree.items_[i].first;
        doc << Circle(Point(point[0], point[1]), 1, Fill(), Stroke(1, Color(200, 250, 150)));

        // circle border
        int id = vptree.items_[i].second;
        const auto& node = vptree.nodes_[id];
        doc << Circle(Point(point[0], point[1]), node.threshold, Fill(), Stroke(1, Color(100, 200, 120)));

        // doc << Circle(Point(point[0], point[1]), 1, Fill(Color(100, 200, 120)), Stroke(1, Color(200, 250, 150)));
    }
}   

void drawVPTree(int node, Polygon area, Color color) 
{
    if(node == vpt::Node::Leaf)
        return;

    var inside = area;
    var outside = area;
    if (node.isRoot) 
    {
        var circle = new vppaper.Path.Circle(node.vp, node.mu);
        node.circle = circle;
        circle.strokeColor = color;
        circle.strokeWidth = strokeWidth;

        trash.push(circle);

        inside = area.intersect(circle);
        outside = area.subtract(circle);
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

int main()
{   
    int dx = 1000, dy = 1000;
    int scale = 1;
    Dimensions dimensions(dx, dy);
    Document doc("vptree.svg", Layout(dimensions, Layout::BottomLeft, scale, Point(0, 0)));

    // // Red image border.
    Polygon border(Fill(Color(255, 255, 255)), Stroke());
    border << Point(0, 0) << Point(dimensions.width, 0) << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    doc << border;

    // auto arc = Arc(Point(dx/2, dy/2), 250, 0, 90, Fill(Color(255, 0, 0)), Stroke(5, Color(0, 255, 0))) ;
    auto arc = Arc(Point(dx/2, dy/2), 250, 0, 270, Fill(), Stroke(5, Color(0, 255, 0))) ;
    // std::cout << arc.toString(Layout(dimensions, Layout::BottomLeft, scale, Point(0, 0))) << std::endl;
    // doc << arc;
    // // Long notation.  Local variable is created, children are added to varaible.
    // LineChart chart(5.0);
    // Polyline polyline_a(Stroke(.5, Color::Blue));
    // Polyline polyline_b(Stroke(.5, Color::Aqua));
    // Polyline polyline_c(Stroke(.5, Color::Fuchsia));
    // polyline_a << Point(0, 0) << Point(10, 30) << Point(20, 40) << Point(30, 45) << Point(40, 44);
    // polyline_b << Point(0, 10) << Point(10, 22) << Point(20, 30) << Point(30, 32) << Point(40, 30);
    // polyline_c << Point(0, 12) << Point(10, 15) << Point(20, 14) << Point(30, 10) << Point(40, 2);
    // chart << polyline_a << polyline_b << polyline_c;
    // doc << chart;

    // // Condensed notation, parenthesis isolate temporaries that are inserted into parents.
    // doc << (LineChart(Dimensions(65, 5))
    //     << (Polyline(Stroke(.5, Color::Blue)) << Point(0, 0) << Point(10, 8) << Point(20, 13))
    //     << (Polyline(Stroke(.5, Color::Orange)) << Point(0, 10) << Point(10, 16) << Point(20, 20))
    //     << (Polyline(Stroke(.5, Color::Cyan)) << Point(0, 5) << Point(10, 13) << Point(20, 16)));

    // doc << Circle(Point(80, 80), 20, Fill(Color(100, 200, 120)), Stroke(1, Color(200, 250, 150)));

    // doc << Text(Point(5, 77), "Simple SVG", Color::Silver, Font(10, "Verdana"));

    // doc << (Polygon(Color(200, 160, 220), Stroke(.5, Color(150, 160, 200))) << Point(20, 70)
    //     << Point(25, 72) << Point(33, 70) << Point(35, 60) << Point(25, 55) << Point(18, 63));

    // doc << Rectangle(Point(70, 55), 20, 15, Color::Yellow);

    int nb = 1000;
    auto points = std::vector<std::vector<double>>(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = drand48() * dx;
        double y = drand48() * dy;
        points[i] = std::vector<double>({x,y,0});
    }
    
    vpt::VpTree t1(points); // create a tree
    
    std::vector<double> distances;
    std::vector<int> indices;
    std::tie(distances, indices) = t1.getNearestNeighbors({ 0, 0, 0 }, 3); // find 3 neighbors closest to the given point
    
    std::cout << distances[0] << "\n"; // prints 0
    std::cout << indices[0] << "\n"; // prints 1
    
    auto batch = t1.getNearestNeighborsBatch({{0, 0, 0}, {1, 1, 1}, {0.5, 0.5, 0.5}}, 3); // split the work between threads
    std::cout << (int)(batch.first[0] == distances) << std::endl; // true
    std::cout << (int)(batch.second[0] == indices) << std::endl; // true

    draw(t1, doc);
    doc.save();

    return 0;
}