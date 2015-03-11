#include "WireCell/CellProperties.h"

WireCell::CellProperties::CellProperties(int id, float area, 
					 const std::vector<float> &center)
    : id(id)
    , area(area)
    , center(center)
{
}
WireCell::CellProperties::~CellProperties()
{
}
