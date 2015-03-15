#include "WireCell/Wire.h"

WireCell::Wire::Wire(int id, 
		     const Point& point1,
		     const Point& point2)
    : id(id)
    , point1(point1)
    , point2(point2)
{
}
WireCell::Wire::~Wire()
{
}

