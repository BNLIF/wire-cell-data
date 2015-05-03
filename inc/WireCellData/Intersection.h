#ifndef WIRECELLDATA_INTERSECTION
#define WIRECELLDATA_INTERSECTION

#include "WireCellData/Vector.h"

namespace WireCell {

    /// Determine if ray hits a square.
    int hit_square(int axis0, 
		   const Vector& bmin, const Vector& bmax, 
		   const Vector& point, const Vector& dir,
		   Vector& hit1, Vector& hit2);

    /// Determine if a ray hits a box.
    int box_intersection(const Vector& bmin, const Vector& bmax, 
			 const Vector& point, const Vector& dir,
			 Vector& hit1, Vector& hit2);

}

#endif
