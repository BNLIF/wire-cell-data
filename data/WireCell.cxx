#include "WireCell/WireCell.h"

WireCell::Id::Id(int cellid, const std::vector<int>& wireids)
    : cell(cellid)
    , wire(wireids)
{
}
WireCell::Id::~Id()
{
}

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
