/** A 3D point of floats
 *
 * See also WCP::Vector.
 *
 * Use Vector unless you need the smaller memory footprint of floats.
 */ 

#ifndef WCPData_Point_h
#define WCPData_Point_h

#include "WCPData/Vector.h"

#include <map>
#include <vector>

namespace WCP {

    typedef D3Vector<float> Point;
    
    typedef std::vector<WCP::Point> PointVector;
    typedef std::pair<WCP::Point, float> PointValue;
    typedef std::vector<WCP::PointValue> PointValueVector;

    /* bool operator==(const WCP::Point& a, const WCP::Point& b){ */
    /*   if (a.x==b.x && a.y==b.y && a.z==b.z){ */
    /* 	return true; */
    /*   }else{ */
    /* 	return false; */
    /*   } */
    /* } */

}
#endif
