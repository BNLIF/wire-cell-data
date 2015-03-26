#include "WireCellData/Wire.h"

#include <algorithm>    // std::sort

WireCellData::Wire::Wire(int ident, 
			 WirePlaneType_t plane,
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


static bool _by_planeindex(const WireCellData::Wire* a, const WireCellData::Wire* b)
{
    if (a->plane < b->plane) return true;
    return a->index < b->index;
}

void WireCellData::sort_by_planeindex(WireCellData::WireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_planeindex);
}

static bool _by_channel(const WireCellData::Wire* a, const WireCellData::Wire* b)
{
    return a->channel < b->channel;
}

void WireCellData::sort_by_channel(WireCellData::WireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_channel);
}
