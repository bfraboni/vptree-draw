#pragma once

#include "simple_svg_1.0.0.hpp" 
#include "cavc/polylinecombine.hpp"

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
      typedef std::array<double, 6>  Edge;
         
      CavcPoly(const cavc::Polyline<double>& poly, Fill const & fill, std::vector<Edge>&_edgeBuffer, Stroke const & stroke = Stroke())
            : Shape(fill, stroke), poly(poly), edgeBuffer(_edgeBuffer) { }
  
       void resetEdgeBuffer() { edgeBuffer.clear(); }
      
        std::string toString(Layout const & layout) const
        {
            std::stringstream ss;

            ss << elemStart("path");
            ss << "d=\"";
            
            for( int i = 0; i < (int)poly.size(); i++)
            {
                const auto& p1 = poly.vertexes()[i];
                const auto& p2 = poly.vertexes()[(i+1)%poly.size()];
                  
                Edge e = { p1.x(), p1.y(), p1.bulge(), p2.x(), p2.y(), p2.bulge()};
                Edge ee = { p2.x(), p2.y(), -p2.bulge(), p1.x(), p1.y(), -p1.bulge()};
                if ((std::find(edgeBuffer.begin(), edgeBuffer.end(), e) == edgeBuffer.end()) && 
                    (std::find(edgeBuffer.begin(), edgeBuffer.end(), ee) == edgeBuffer.end()))
                {
                    if(!p1.bulgeIsZero())
                    {
                        // arc if bulge != 0
                        auto arc = arcRadiusAndCenter(p1, p2);
                        double radius = arc.radius;
                        Point ps(p1.pos().x(), p1.pos().y());
                        Point pe(p2.pos().x(), p2.pos().y());

                        ss << "M";
                        ss << translateX(ps.x, layout) << "," << translateY(ps.y, layout);
                        ss << "A";
                        ss << translateScale(radius, layout) << "," << translateScale(radius, layout);
                        ss << " 0 ";
                        // always small arcs -> large arc flag to 0
                        ss << "0,";                           
                        // sweep is 0 if bulge is positive, 1 if negative
                        ss << (p1.bulge() > 0 ? "0 " : "1 "); 
                        ss << translateX(pe.x, layout) << "," << translateY(pe.y, layout);
                    }
                    else
                    {
                        // line if bulge == 0
                        ss << "M";
                        ss << translateX(p1.pos().x(), layout) << "," << translateY(p1.pos().y(), layout) << " ";
                        ss << translateX(p2.pos().x(), layout) << "," << translateY(p2.pos().y(), layout) << " ";
                    }
                
                    edgeBuffer.push_back(e);
                }
            }
            ss << "\" ";
            ss << fill.toString(layout) << stroke.toString(layout) << emptyElemEnd();

            return ss.str();
        }
        void offset(Point const & offset) { }
    private:
        const cavc::Polyline<double>& poly;
        std::vector<Edge> &edgeBuffer;
    };
}