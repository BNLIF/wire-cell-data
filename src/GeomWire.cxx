#include "WireCellData/GeomWire.h"

#include <algorithm>    // std::sort


using namespace WireCell;

GeomWire::GeomWire(int ident,
		   WirePlaneType_t plane,
		   int index,
		   int channel,
		   const Point& point1,
		   const Point& point2,
		   char segment,
		   char face,
		   short apa)
    : _ident(ident)
    , _plane(plane)
    , _index(index)
    , _channel(channel)
    , _point1(point1)
    , _point2(point2)
    , _segment(segment)
    , _face(face)
    , _apa(apa)
{
}

GeomWire::~GeomWire()
{
}

std::ostream & WireCell::operator<<(std::ostream &os, const GeomWire& gw)
{
    return os << "<WireCell::GeomWire "
	      << "id:" << gw.ident() << " "
	      << "ch:" << gw.channel() << " "
	      << gw.apa() << "/"
	      << gw.face() << "/"
	      << gw.plane() << "/"
	      << gw.segment() << ">";
}


bool _by_planeindex(const GeomWire* a, const GeomWire* b)
{
    if (a->plane() < b->plane()) return true;
    return a->index() < b->index();
}

void WireCell::sort_by_planeindex(GeomWireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_planeindex);
}

static bool _by_channel(const GeomWire* a, const GeomWire* b)
{
    return a->channel() < b->channel();
}

void WireCell::sort_by_channel(GeomWireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_channel);
}

