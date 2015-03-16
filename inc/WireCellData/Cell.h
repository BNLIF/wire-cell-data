#ifndef WireCellData_Cell_h
#define WireCellData_Cell_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

namespace WireCellData {

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
