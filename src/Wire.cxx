#include "WireCellData/Wire.h"

WireCellData::Wire::Wire(int ident, 
			 int plane, 
			 int index,
			 int channel,
			 const Point& point1,
			 const Point& point2)
    : ident(ident)
    , plane(plane)
    , index(index)
    , channel(channel)
    , point1(point1)
    , point2(point2)
{
}
WireCellData::Wire::~Wire()
{
}

