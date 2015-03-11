#include "WireCell/Id.h"

WireCell::Id::Id(int cellid, const std::vector<int>& wireids)
    : cell(cellid)
    , wire(wireids)
{
}
WireCell::Id::~Id()
{
}

