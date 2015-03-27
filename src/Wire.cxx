#include "WireCellData/Wire.h"

#include <algorithm>    // std::sort

using namespace WireCell;

Wire::Wire(int ident, 
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
Wire::~Wire()
{
}


static bool _by_planeindex(const Wire* a, const Wire* b)
{
    if (a->plane < b->plane) return true;
    return a->index < b->index;
}

void WireCell::sort_by_planeindex(WireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_planeindex);
}

static bool _by_channel(const Wire* a, const Wire* b)
{
    return a->channel < b->channel;
}

void WireCell::sort_by_channel(WireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_channel);
}
