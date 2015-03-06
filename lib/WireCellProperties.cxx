#include "WireCellProperties.h"

WireCell::Properties::Properties(int id, float area_cm2, 
				 const std::vector<float> &center_cm)
    : cell(id)
    , area_cm2(area_cm2)
    , center_cm(center_cm)
{
}
WireCell::Properties::~Properties()
{
}
