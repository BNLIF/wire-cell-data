#include "WireCellData/Intersection.h"

#include <iostream>
using namespace std;

using namespace WireCell;



/// Return 0 if no hit, 1 if hit1, 2 if hit2, 3 if both
int WireCell::hit_square(int axis0, 
			 const D3FloatVector& bmin, 
			 const D3FloatVector& bmax, 
			 const D3FloatVector& point, 
			 const D3FloatVector& dir,
			 D3FloatVector& hit1, 
			 D3FloatVector& hit2)
{
    int hitmask = 0;

    if (0 == dir[axis0]) {
	return hitmask;
    }

    int axis1 = (axis0 + 1)%3;
    int axis2 = (axis1 + 1)%3;

    { // toward the min intercept
	float intercept = bmin[axis0];
	float scale = (intercept - point[axis0])/dir[axis0];

	float one = point[axis1] + scale*dir[axis1];
	float two = point[axis2] + scale*dir[axis2];

	cerr << "MIN: " << axis0 << " scale=" << scale
	     << " one:" << one << " in:[" << bmin[axis1] << "," << bmax[axis1] << "]"
	     << " two:" << two << " in:[" << bmin[axis2] << "," << bmax[axis2] << "]"
	     <<endl;
    
	if (bmin[axis1] <= one && one <= bmax[axis1] &&
	    bmin[axis2] <= two && two <= bmax[axis2]) { 
	    hitmask |= 1;
	    hit1[axis0] = intercept;
	    hit1[axis1] = one;
	    hit1[axis2] = two;
	    cerr << "HIT1:" << hit1 << endl;
	}
    }

    { // toward the max intercept
	float intercept = bmax[axis0];
	float scale = (intercept - point[axis0])/dir[axis0];

	float one = point[axis1] + scale*dir[axis1];
	float two = point[axis2] + scale*dir[axis2];

	cerr << "MAX: " << axis0 << " scale=" << scale
	     << " one:" << one << " in:[" << bmin[axis1] << "," << bmax[axis1] << "]"
	     << " two:" << two << " in:[" << bmin[axis2] << "," << bmax[axis2] << "]"
	     <<endl;

	if (bmin[axis1] <= one && one <= bmax[axis1] && 
	    bmin[axis2] <= two && two <= bmax[axis2]) {
	    hitmask |= 2;
	    hit2[axis0] = intercept;
	    hit2[axis1] = one;
	    hit2[axis2] = two;
	    cerr << "HIT2:" << hit2 << endl;
	}
    }
		
    return hitmask;
}
		
/// Return 0 if no hit, 1 if hit1, 2 if hit2, 3 if both
int WireCell::box_intersection(const D3FloatVector& bmin, 
			       const D3FloatVector& bmax, 
			       const D3FloatVector& point, 
			       const D3FloatVector& dir,
			       D3FloatVector& hit1, 
			       D3FloatVector& hit2)
{
    std::vector<D3FloatVector> results;

    for (int axis=0; axis<3; ++axis) {
	D3FloatVector h1, h2;
	int got = hit_square(axis, bmin, bmax, point, dir, h1, h2);
	if (got&1) {
	    results.push_back(h1);
	}
	if (got&2) {
	    results.push_back(h2);
	}
    }

    if (results.size() > 2) {
	cerr << "ERROR: crossed box more than twice" << endl;
    }

    int hitmask = 0;
    for (int ind=0; ind < results.size() && ind < 2; ++ind) {
	D3FloatVector& hit = results[ind];
	D3FloatVector hitdir = hit - point;
	float dot = hitdir.norm().dot(dir);

	cerr << "DOT: " << dot << " " << hitdir<< dir << hit << point << endl;

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

