#ifndef WireCellData_Point_h
#define WireCellData_Point_h

#include <map>
#include <vector>

namespace WireCell {
    struct Point {
	Point(float x=0, float y=0, float z=0) : x(x), y(y), z(z) { }
	float x, y, z;
    };
    typedef std::vector<WireCell::Point> PointVector;

    struct PointC {
    PointC(float x=0, float y=0, float z=0, float charge=0) : x(x), y(y), z(z), charge(charge){ }
      float x, y, z, charge;
    };
    typedef std::vector<WireCell::PointC> PointCVector;
}
#endif
