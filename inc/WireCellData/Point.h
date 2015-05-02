#ifndef WireCellData_Point_h
#define WireCellData_Point_h

#include "WireCellData/Vector.h"

#include <map>
#include <vector>

namespace WireCell {

    typedef D3Vector<float> Point;

    typedef std::vector<WireCell::Point> PointVector;
    typedef std::pair<WireCell::Point, float> PointValue;
    typedef std::vector<WireCell::PointValue> PointValueVector;

}
#endif
