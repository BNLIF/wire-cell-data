#include "WCPData/Vector.h"
#include "WCPData/Intersection.h"

#include <iostream>
#include <cmath>
#include <map>

using namespace WCP;
using namespace std;


int main()
{
    D3Vector<int> a( 3 , 4 , 5 ) , b ( 4 , 3 , 5 ) , c( -5 , -12 , -13 ) ;
    cout << "a . b : " << a.dot( b ) << endl ;
    cout << "a x b : " << a.cross( b ) << endl ;
    cout << "a . b x c : " << a.triplescal( b , c ) << endl ;
    cout << "a x b x c : " << a.triplevec( b , c ) << endl ;

    Vector b1(0,0,0), b2(1,1,1), ray(1,1,1);
    ray = ray.norm();

    for (double x = -1.1; x <= 1; x+=0.5) {
	for (double y = -1.1; y <= 1; y+=0.5) {
	    for (double z = -1.0; z <= 1; z+=0.5) {
		Vector point(x,y,z);
		Vector hit1(-111,-111,-111), hit2(-222,-222,-222);

		int hitmask = box_intersection(b1, b2, point, ray, hit1, hit2);
		cerr << "RESULT: " 
		     << hitmask << " "
		     << point << " -> " << ray 
		     << " hit1=" << hit1 << " hit2=" << hit2 << endl;
	    }
	}
    }

    //pair<float,float> pp = box_interesect(b1, b2, off, pitch);
    //cout << b1 << " -> " << b2 << " from " << off << " to " << pitch << endl;
    //cout << "distances: " << pp.first << ", " << pp.second << endl;


    return 0 ;
}
