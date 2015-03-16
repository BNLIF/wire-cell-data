#include "WireCellData/Cell.h"

WireCellData::Cell::Cell(int ident, float area, const Point& center, 
			 const PointVector& boundary)
    : ident(ident)
    , area(area)
    , center(center)
    , boundary(boundary)
{
}
WireCellData::Cell::~Cell()
{
}
