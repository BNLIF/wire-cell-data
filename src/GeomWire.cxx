#include "WCPData/GeomWire.h"

#include <algorithm>    // std::sort


using namespace WCP;

/*
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
*/

GeomWire::GeomWire(unsigned int ident,
		   WirePlaneType_t plane,
		   int index,
		   int channel,
		   const Point& point1,
		   const Point& point2,
		   char segment,
		   int face,
		   short apa,
		   short cryo)
    : _ident(ident)
    , _plane(plane)
    , _index(index)
    , _channel(channel)
    , _point1(point1)
    , _point2(point2)
    , _segment(segment)
    , _face(face)
    , _cryo(cryo)
    , _apa(apa)
{
}


GeomWire::~GeomWire()
{
}

std::ostream & WCP::operator<<(std::ostream &os, const GeomWire& gw)
{
    char plane_name[] = {'U', 'V', 'W'};
    return os << "<WCP::GeomWire"
	      << " " << plane_name[gw.iplane()] << " "
	      << "id:" << gw.ident() << " "
	      << "ch:" << gw.channel() << " "
	      << "f:" << gw.face() << " "
	      << "i:" << gw.index() << " "
	      << "s:" << gw.segment() << " "
	      << "a:" << gw.apa() << " "
	      << "c:" << gw.cryo() << ">";
}


bool _by_planeindex(const GeomWire* a, const GeomWire* b)
{
    if (a->plane() < b->plane()) return true;
    return a->index() < b->index();
}

void WCP::sort_by_planeindex(GeomWireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_planeindex);
}

bool _by_ident(const GeomWire* a, const GeomWire* b)
{
    if (a->ident() < b->ident()) return true;
    return a->ident() < b->ident();
}

void WCP::sort_by_ident(GeomWireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_ident);
}

static bool _by_channel(const GeomWire* a, const GeomWire* b)
{
    return a->channel() < b->channel();
}

void WCP::sort_by_channel(GeomWireSelection& ws)
{
    std::sort(ws.begin(), ws.end(), _by_channel);
}

double WCP::GeomWire::length() const {
  double lsquare = TMath::Power((_point1.x - _point2.x),2) 
    + TMath::Power((_point1.y - _point2.y),2)
    + TMath::Power((_point1.z - _point2.z),2);
  return TMath::Sqrt(lsquare);
}
