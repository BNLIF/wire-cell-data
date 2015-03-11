#include "CellProperties.h"

WireCell::CellProperties::CellProperties(int id, float area_cm2, 
					 const std::vector<float> &center_cm)
    : id(id)
    , area_cm2(area_cm2)
    , center_cm(center_cm)
{
}
WireCell::CellProperties::~CellProperties()
{
}
