#ifndef WireCellData_Wire_h
#define WireCellData_Wire_h

#include <vector>

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

namespace WireCellData {

/** 
    Collect information about one wire.
 */
struct Wire {

    Wire(int ident = 0,  
	 int plane = 0,  
	 int index = -1,
	 int channel = -1,
	 const Point& point1 = Point(),
	 const Point& point2 = Point());
    ~Wire();

    /// globally unique ID (zero is illegal)
    int ident;
    // the plane number (zero is illegal)
    int plane;
    // index into ordered sequence of wires in a plane
    int index;
    // electronics channel connected
    int channel;
    // End points of the wire
    Point point1, point2;
};

    typedef std::vector<Wire> WireCollection;

/**
   A wire triple is simply one wire from each plane identified by their ID numbers.
 */
struct WireTriple {
    int u, v, z;
    WireTriple(int u, int v, int z) : u(u), v(v), z(z) {}
};


}
#endif
