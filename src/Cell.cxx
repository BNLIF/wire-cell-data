#include "WireCellData/Cell.h"

using namespace WireCell;

Cell::Cell(int ident, float area, const Point& center, 
		     const PointVector& boundary)
    : ident(ident)
    , area(area)
    , center(center)
    , boundary(boundary)
{
}
Cell::~Cell()
{
}
