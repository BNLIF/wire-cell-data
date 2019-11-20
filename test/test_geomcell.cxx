#include "WCPData/GeomCell.h"

#include <iostream>

using namespace std;
using namespace WCP;

int main() {

    PointVector pv1;
    pv1.push_back(Point(0,+1,+1));
    pv1.push_back(Point(0,-1,+1));
    pv1.push_back(Point(0,-1,-1));
    pv1.push_back(Point(0,+1,-1));

    GeomCell gc1(0, pv1);

    Point center = gc1.center();
    if (center.x != 0.0 ||
	center.y != 0.0 ||
	center.z != 0.0) {

	cerr << "Wrong center from GeomCell"
	     << " " << center.x
	     << " " << center.y
	     << " " << center.z
	     << endl;
	exit(1);
    }


    float area = gc1.cross_section();
    if (area != 4.0) {
	cerr << "Wrong area from GeomCell " << area << endl;
	exit(1);
    }

    return 0;
}
