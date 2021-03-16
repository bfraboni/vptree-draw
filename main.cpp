#include <iostream> 
#include <random>
#include <set>
#include <cassert>
#include <set>

#include "vptree.hpp" 
#include "simple_svg_1.0.0.hpp" 
#include "cavc/polylinecombine.hpp"

using namespace svg;

struct Cell
{
    std::vector<cavc::Polyline<double>> shape;
    std::vector<cavc::Polyline<double>> holes;
};

void drawVPTree(
    const vpt::VpTree& vptree, 
    int depth, 
    int node, 
    const Cell& cell,
    std::vector<CavcPoly::Edge> &edgeBuffer, 
    svg::Document& doc
)
{
    if(node == vpt::VpTree::Node::Leaf)
        return;

    // get cell vantage circle info
    const auto& n = vptree.nodes_[node];
    double radius = n.threshold;
    double cx = vptree.items_[n.item].first[0];
    double cy = vptree.items_[n.item].first[1];
    
    // hue rotated color 
    svg::Color color = rotate(svg::Color(255,120,80), depth*30);
    std::cout << "depth: " << depth << " " << node << " " << n.left << " " << n.right << std::endl;
    std::cout << "circle: " << cx << " " << cy << " " << radius << std::endl;

    // on a leaf
    if( radius <= 0 )
    {
        // draw vantage point
        doc << svg::Circle(svg::Point(cx, cy), 4, Fill(color), Stroke());
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
        for(int i = 0; i < (int)cell.shape.size(); ++i)
        {
            cavc::CombineResult<double> inter = combinePolylines(cell.shape[i], circle, cavc::PlineCombineMode::Intersect);
            std::cout << "intersect: " << inter.remaining.size() << " remains " << inter.subtracted.size() << " holes " << std::endl;
            if( inter.remaining.size() > 0 )
            {
                inside.shape.insert(inside.shape.end(), inter.remaining.begin(), inter.remaining.end());
            }
        }

        // compute exclusion with enclosing path
        Cell outside;
        for(int i = 0; i < (int)cell.shape.size(); ++i)
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
                outside.holes.insert(outside.holes.end(), exclu.subtracted.begin(), exclu.subtracted.end());
            }
        }

        // draw circle center
        doc << svg::Circle(svg::Point(cx, cy), 4, Fill(color), Stroke());

        // draw left subtree
        drawVPTree(vptree, depth+1, n.left, inside, edgeBuffer, doc);

        // draw right subtree
        drawVPTree(vptree, depth+1, n.right, outside, edgeBuffer, doc);

        // draw left cell (inside)
        for(int i = 0; i < (int)inside.shape.size(); ++i)
        {
            doc << svg::CavcPoly(inside.shape[i], Fill(), edgeBuffer, Stroke(2, color));
        }

        // draw right cell (outside)
        for(int i = 0; i < (int)outside.shape.size(); ++i)
        {
            doc << svg::CavcPoly(outside.shape[i], Fill(), edgeBuffer, Stroke(2, color));
        }    
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

void initbunny(int dx, int dy, std::vector<std::vector<double>>& points)
{
    std::vector<double> data;
    std::ifstream file("bunny.dat");
    double x, y;
    while(file)
    {
        file >> x >> y;
        data.push_back(x);
        data.push_back(y);
    }
    if(data.empty()) exit(1);
    int nb = data.size()/2;
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        points[i] = std::vector<double>({data[2*i]*dx,data[2*i+1]*dy,0});
    }
}

void initrandom(int dx, int dy, std::vector<std::vector<double>>& points)
{
    int nb = 1000;
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = rand1D() * dx;
        double y = rand1D() * dy;
        points[i] = std::vector<double>({x,y,0});
    }
}

void initdiag(int dx, int dy, std::vector<std::vector<double>>& points)
{
    int nb = 1000;
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = rand1D() * dx;
        points[i] = std::vector<double>({x,x,0});
    }
}

static double gaussian(double x, double m, double s)
{
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * std::exp(-0.5 * a * a);
}

void initgauss(int dx, int dy, std::vector<std::vector<double>>& points)
{
    int nb = 1000;
    points.resize(nb);
    for(int i = 0; i < nb; ++i)
    {
        double x = rand1D();
        double y = gaussian(3.0*x, 1.5, 0.4);
        points[i] = std::vector<double>({x*dx,y*dy,0});
    }
}

void initsquare(int dx, int dy, std::vector<std::vector<double>>& points)
{
    int nb = 1000;
    points.resize(nb);
    for(int i = 0; i < nb/2; ++i)
    {
        double x = rand1D();
        points[2*i] = std::vector<double>({x*dx,x*x*dy,0});
        points[2*i+1] = std::vector<double>({x*dx,(1.f-x*x)*dy,0});
    }
}

int main()
{   
    int dx = 1000, dy = 1000;
    int scale = 1;
    Dimensions dimensions(dx, dy);
    Layout layout(dimensions, Layout::BottomLeft, scale, Point(0, 0));
    Document doc("vpbunny.svg", layout);

    // background white rectangle
    Polygon border(Fill(Color(255, 255, 255)), Stroke());
    border << Point(0, 0) << Point(dimensions.width, 0) << Point(dimensions.width, dimensions.height) << Point(0, dimensions.height);
    doc << border;

    // init random points
    std::vector<std::vector<double>> points;
    initbunny(dx, dy, points);
    // initgauss(dx, dy, points);
    // initdiag(dx, dy, points);
    // initsquare(dx, dy, points);

    // create a vantage points tree
    vpt::VpTree vptree(points);

    // init cell representing the canvas
    cavc::Polyline<double> shape1;
    shape1.addVertex(1, 1, 0); 
    shape1.addVertex(dx-1, 1, 0); 
    shape1.addVertex(dx-1, dy-1, 0); 
    shape1.addVertex(1, dy-1, 0);
    shape1.isClosed() = true;
    
    Cell cell;
    cell.shape.push_back(shape1);
    
    // init edge structure
    std::vector<CavcPoly::Edge> buffer;
  
    // draw VP tree
    drawVPTree(vptree, 0, 0, cell, buffer, doc);
    
    // save svg
    doc.save();

    return 0;
}
