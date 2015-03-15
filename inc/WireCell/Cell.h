#ifndef Cell_h
#define Cell_h

#include "WireCell/Units.h"
#include "WireCell/Point.h"

namespace WireCell {

struct Cell {

    Cell(int ident = 0, float area=0.0 * units::centimeter2, 
	 const Point& center=Point(), 
	 const PointVector& boundary = PointVector());

    ~Cell();

    int ident;
    float area;
    Point center;
    PointVector boundary;
};
}
#endif
