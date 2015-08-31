#include "WireCellData/Plane.h"

using namespace WireCell;

Plane::Plane(Point &p1, Point &p2, Point &p3)
  : p1(p1)
  , p2(p2)
  , p3(p3)
{
}

Plane::~Plane(){
}
