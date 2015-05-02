#ifndef WIRECELLDATA_INTERSECTION
#define WIRECELLDATA_INTERSECTION

#include "WireCellData/Vector.h"

namespace WireCell {

    /// Determine if ray hits a square.
    int hit_square(int axis0, 
		   const D3FloatVector& bmin, 
		   const D3FloatVector& bmax, 
		   const D3FloatVector& point, 
		   const D3FloatVector& dir,
		   D3FloatVector& hit1, 
		   D3FloatVector& hit2);

    /// Determine if a ray hits a box.
    int box_intersection(const D3FloatVector& bmin, 
			 const D3FloatVector& bmax, 
			 const D3FloatVector& point, 
			 const D3FloatVector& dir,
			 D3FloatVector& hit1, 
			 D3FloatVector& hit2);

}

#endif
