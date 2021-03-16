#pragma once

#include <cmath>
#include <limits>

#define EPSILON 1e-6f

namespace geo
{
    struct vec2
    {
        float x, y;
        float& operator[](int i) {return *(&x+i);}
        float operator[](int i) const {return *(&x+i);}
        float& operator()(int i) {return *(&x+i);}
        float operator()(int i) const {return *(&x+i);}
        
        friend std::ostream& operator<< (std::ostream& out, const vec2& v)
        {
            return out << "(" << v.x << ";" << v.y << ")";
        }

        vec2() : x(), y() {}
        vec2(float v) : x(v), y(v) {}
        vec2(float x, float y) : x(x), y(y) {}
        vec2(const vec2& a, const vec2& b) : x(b.x-a.x), y(b.y-a.y) {}

        bool operator==(const vec2& other) const { return x == other.x && y == other.y; }
        bool operator!=(const vec2& other) const { return !(*this==other); }

        vec2 operator+(const vec2& other) const { return vec2(x+other.x, y+other.y); }
        vec2 operator-(const vec2& other) const { return vec2(x-other.x, y-other.y); }
        vec2 operator*(const vec2& other) const { return vec2(x*other.x, y*other.y); }
        vec2 operator/(const vec2& other) const { return vec2(x/other.x, y/other.y); }

        vec2 operator+(float k) const { return vec2(x+k, y+k); }
        vec2 operator-(float k) const { return vec2(x-k, y-k); }
        vec2 operator*(float k) const { return vec2(x*k, y*k); }
        vec2 operator/(float k) const { return vec2(x/k, y/k); }

        vec2& operator+=(const vec2& other) { x+=other.x; y+=other.y; return *this; }
        vec2& operator-=(const vec2& other) { x-=other.x; y-=other.y; return *this; }
        vec2& operator*=(const vec2& other) { x*=other.x; y*=other.y; return *this; }
        vec2& operator/=(const vec2& other) { x/=other.x; y/=other.y; return *this; }

        vec2& operator+=(float k) { x+=k; y+=k; return *this; }
        vec2& operator-=(float k) { x-=k; y-=k; return *this; }
        vec2& operator*=(float k) { x*=k; y*=k; return *this; }
        vec2& operator/=(float k) { x/=k; y/=k; return *this; }

        static vec2 ones() { return vec2(1,1); }
        static vec2 zero() { return vec2(0,0); }
    };

    vec2 operator*(float k, const vec2& v) { return vec2(v.x*k, v.y*k); }

    vec2 lerp(const vec2& a, const vec2& b, const float t) 
    {
        if( t == 1 ) return b;
        if( t == 0 ) return a;
        return a + t * (b - a);
    }

    float distance(const vec2& a, const vec2& b)
    {
        return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
    } 

    vec2 min( const vec2& a, const vec2& b ) 
    { 
        return vec2( std::min(a.x, b.x), std::min(a.y, b.y) ); 
    }

    vec2 max( const vec2& a, const vec2& b ) 
    { 
        return vec2( std::max(a.x, b.x), std::max(a.y, b.y) ); 
    }

    typedef vec2 Point;
    typedef vec2 Vector;

    float dot(const Vector& a, const Vector& b) { return a.x * b.x + a.y * b.y; }
    float length2(const Vector& a) { return dot(a, a); }
    float length(const Vector& a) { return std::sqrt(length2(a)); }
    Vector normalize(const Vector& a) { const float l = length(a); return a / l; }

    struct Box
    {
        Point pmin = Point(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        Point pmax = Point(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
        Box(){}
        Box( const Point& p) : pmin(p), pmax(p) {}
        Box( const Point& pmin, const Point& pmax) : pmin(pmin), pmax(pmax) {}
        Box( const Box& a, const Box& b ) : pmin(min(a.pmin, b.pmin)), pmax(max(a.pmax, b.pmax)) {}
        bool empty() const { return ((pmax.x <= pmin.x) && (pmax.y <= pmin.y)); }
        void insert( const Point& p ) { pmin = min(pmin, p); pmax = max(pmax, p); }
        void insert( const Box& b ) { pmin = min(pmin, b.pmin); pmax = max(pmax, b.pmax); }
        Point center() const {return pmin + (pmax - pmin) / 2.f;}
        bool contains(const Point& p) const 
        { 
            return pmin.x <= p.x && p.x < pmax.x && pmin.y <= p.y && p.y < pmax.y; 
        }
    };

    struct Sphere 
    {
        Point pos = Point(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        float r = -1.f;
        Sphere(){}
        Sphere( const Point& p, const float r = EPSILON) : pos(p), r(r) {}
        Sphere( const Sphere& a, const Sphere& b) 
        {
            if( a.pos == b.pos )
            {
                pos = a.pos;
                r = std::max(a.r, b.r);
            }
            else
            {
                // nouveau centre = point milieu des deux extremités
                Vector l = normalize(Vector(b.pos, a.pos));
                Vector p1 = b.pos - l * (b.r + EPSILON);
                Vector p2 = a.pos + l * (a.r + EPSILON);
                pos = (p1 + p2) * 0.5f;
                // nouveau rayon 
                r = length(p2 - p1) * 0.5f; 
            }
        }
        bool empty() const { return r <= 0; }
        bool contains( const Point& p ) const { return distance(p, pos) <= r; }
        Point center() const { return pos; }
        bool operator==(const Sphere& s) const { return pos == s.pos && r == s.r; }
        void insert( const Point& p ) { insert(Sphere(p, 0)); }
        void insert( const Sphere& s )
        {
            if(s.empty() || *this == s) return;

            if(empty())
            {
                pos = s.pos;
                r = s.r;
                return;
            }
            
            // nouveau centre = point milieu des deux extremités
            Vector l = normalize(Vector(s.pos, pos));
            Vector p1 = s.pos - l * (s.r + EPSILON);
            Vector p2 = pos + l * (r + EPSILON);
            pos = (p1 + p2) * 0.5f;
            // nouveau rayon 
            r = length(p2 - p1) * 0.5f; 
        }
    };
}