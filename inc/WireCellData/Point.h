/** A 3D point of floats
 *
 * See also WireCell::Vector.
 *
 * Use Vector unless you need the smaller memory footprint of floats.
 */ 

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

    /* bool operator==(const WireCell::Point& a, const WireCell::Point& b){ */
    /*   if (a.x==b.x && a.y==b.y && a.z==b.z){ */
    /* 	return true; */
    /*   }else{ */
    /* 	return false; */
    /*   } */
    /* } */

}
#endif
