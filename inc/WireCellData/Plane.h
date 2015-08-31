#ifndef WireCellData_Plane_h
#define WireCellData_Plane_h

#include "WireCellData/Line.h"

namespace WireCell {
  class Plane{
  public:
    Plane(Point &p1, Point &p2, Point &p3);
    ~Plane();
  protected:
    Point p1;
    Point p2;
    Point p3;
  };
}

#endif
