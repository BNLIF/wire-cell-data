#include "WireCell/Wire.h"

WireCell::Wire::Wire(int id, float angle,
		     const Point& location)
    : id(id)
    , angle(angle)
    , location(location)
{
}
WireCell::Wire::~Wire()
{
}

