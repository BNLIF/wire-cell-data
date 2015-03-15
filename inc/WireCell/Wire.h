#ifndef WireCell_Wire_h
#define WireCell_Wire_h

#include "WireCell/Units.h"
#include "WireCell/Point.h"

namespace WireCell {

/** 
    Collect information about one wire.
 */
struct Wire {

    Wire(int id = 0, 
	 const Point& point1 = Point(),
	 const Point& point2 = Point());
    ~Wire();

    int id;
    Point point1, point2;
};

/**
   A wire triple is simply one wire from each plane identified by their ID numbers.
 */
struct WireTriple {
  int u, v, y;
  WireTriple(int u, int v, int y) : u(u), v(v), y(y) {}
};


}
#endif
