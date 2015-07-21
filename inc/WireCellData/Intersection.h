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

    /** Return the distance from a point a ray to the intersection of the ray and a plane.
     * 
     * \param is a point on the ray.
     *
     * \param is a vector in the direction of the ray, distance is returned in units of its magnitude.
     *
     * \param is a plane defined as the perpendicular distance from the origin to a point on the plane.
     *
     * Caveat: assumes ray actually hits plane.
     */
    double dist_to_plane(const Vector& point, const Vector& dir, const Vector& plane);

    /** Return the dot product of the directions of two vectors. */
    double directional_dot(const Vector& v1, const Vector& v2);
    
}

#endif
