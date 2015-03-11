#ifndef WireCell_Wire_h
#define WireCell_Wire_h

#include "WireCell/Units.h"
#include "WireCell/Point.h"

namespace WireCell {

struct Wire {

    Wire(int id = 0, float angle = 0.0 * units::degree,
	 const Point& location = Point());
    ~Wire();

    int id;
    float angle;
    Point location;
};
}
#endif
