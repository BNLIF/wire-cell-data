#include "WCPData/Intersection.h"

#include <vector>
#include <iostream>
using namespace std;

using namespace WCP;


double WCP::directional_dot(const Vector& dir1, const Vector& dir2)
{
    return dir1.norm().dot(dir2.norm());
}

double WCP::dist_to_plane(const Vector& point, const Vector& dir, const Vector& plane)
{
    const Vector norm = plane.norm();
    return norm.dot(point - plane) / directional_dot(plane, dir);
}





/// Return 0 if no hit, 1 if hit1, 2 if hit2, 3 if both
int WCP::hit_square(int axis0, 
			 const Vector& bmin, 
			 const Vector& bmax, 
			 const Vector& point, 
			 const Vector& dir,
			 Vector& hit1, 
			 Vector& hit2)
{
    int hitmask = 0;

    if (0 == dir[axis0]) {
	return hitmask;
    }

    int axis1 = (axis0 + 1)%3;
    int axis2 = (axis1 + 1)%3;

    { // toward the min intercept
	double intercept = bmin[axis0];
	double scale = (intercept - point[axis0])/dir[axis0];

	double one = point[axis1] + scale*dir[axis1];
	double two = point[axis2] + scale*dir[axis2];

	// cerr << "MIN: " << axis0 << " scale=" << scale
	//      << " one:" << one << " in:[" << bmin[axis1] << "," << bmax[axis1] << "]"
	//      << " two:" << two << " in:[" << bmin[axis2] << "," << bmax[axis2] << "]"
	//      <<endl;
    
	if (bmin[axis1] <= one && one <= bmax[axis1] &&
	    bmin[axis2] <= two && two <= bmax[axis2]) { 
	    hitmask |= 1;
	    hit1[axis0] = intercept;
	    hit1[axis1] = one;
	    hit1[axis2] = two;
	    //cerr << "HIT1:" << hit1 << endl;
	}
    }

    { // toward the max intercept
	double intercept = bmax[axis0];
	double scale = (intercept - point[axis0])/dir[axis0];

	double one = point[axis1] + scale*dir[axis1];
	double two = point[axis2] + scale*dir[axis2];

	// cerr << "MAX: " << axis0 << " scale=" << scale
	//      << " one:" << one << " in:[" << bmin[axis1] << "," << bmax[axis1] << "]"
	//      << " two:" << two << " in:[" << bmin[axis2] << "," << bmax[axis2] << "]"
	//      <<endl;

	if (bmin[axis1] <= one && one <= bmax[axis1] && 
	    bmin[axis2] <= two && two <= bmax[axis2]) {
	    hitmask |= 2;
	    hit2[axis0] = intercept;
	    hit2[axis1] = one;
	    hit2[axis2] = two;
	    //cerr << "HIT2:" << hit2 << endl;
	}
    }
		
    return hitmask;
}
		
/// Return 0 if no hit, 1 if hit1, 2 if hit2, 3 if both
int WCP::box_intersection(const Vector& bmin, 
			       const Vector& bmax, 
			       const Vector& point, 
			       const Vector& dir,
			       Vector& hit1, 
			       Vector& hit2)
{
    std::vector<Vector> results;

    for (int axis=0; axis<3; ++axis) {
	Vector h1, h2;
	int got = hit_square(axis, bmin, bmax, point, dir, h1, h2);
	if (got&1) {
	    results.push_back(h1);
	}
	if (got&2) {
	    results.push_back(h2);
	}
    }

    if (results.size() > 2) {
	cerr << "ERROR: crossed box " << results.size() << " times." << endl;
    }

    int hitmask = 0;
    for (int ind=0; ind < results.size() && ind < 2; ++ind) {
	Vector& hit = results[ind];
	Vector hitdir = hit - point;
	double dot = hitdir.norm().dot(dir);

	//cerr << "DOT: " << dot << " " << hitdir<< dir << hit << point << endl;

	if (dot > 0) {		// really should be +/- 1 w/in tollerance
	    hit1 = hit;
	    hitmask |= 1;
	}
	else {
	    hit2 = hit;
	    hitmask |= 2;
	}
    }

    return hitmask;
}

