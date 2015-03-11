#include "WireCell/Cell.h"

WireCell::Cell::Cell(int id, float area, const Point& center, 
		     const PointVector& boundary)
    : id(id)
    , area(area)
    , center(center)
    , boundary(boundary)
{
}
WireCell::Cell::~Cell()
{
}
