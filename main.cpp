#include <iostream> 
#include "vptree/vptree.hpp" 
#include "simple_svg_1.0.0.hpp" 

using namespace svg;

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

int main()
{   
    int dx = 1000, dy = 1000;
    int scale = 1;
    Dimensions dimensions(dx, dy);
    Document doc("vptree.svg", Layout(dimensions, Layout::BottomLeft, scale, Point(dx/2/scale, dy/2/scale)));

    // // Red image border.
    // Polygon border(Stroke(1, Color::Red));
    // border << Point(0, 0) << Point(dimensions.width, 0) << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    // doc << border;

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

    auto points = std::vector<std::vector<double>>{
        {0, 0, 1},
        {1, 1, 1},
        {2, 0, 0},
        {-1, -1, 0},
        {10, 0, 5}
    };
    
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