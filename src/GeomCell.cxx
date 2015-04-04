#include "WireCellData/GeomCell.h"

#include <vector>

using namespace std;
using namespace WireCell;

GeomCell::GeomCell(int ident, const PointVector& boundary)
    : _ident(ident)
    , _boundary(boundary)
{
}
GeomCell::~GeomCell()
{
}

std::ostream & WireCell::operator<<(std::ostream &os, const GeomCell& gc)
{
    return os << "<WireCell::GeomCell " << gc.ident() << ">";
}

double GeomCell::cross_section() const
{
    double area = 0.0;
    const size_t npoints = _boundary.size();
    int prev = npoints - 1;

    for (int ind = 0; ind < npoints; ++ind) {
	double z = _boundary.at(prev).z + _boundary.at(ind).z;
	double y = _boundary.at(prev).y - _boundary.at(ind).y;
	area += y*z;
	prev = ind;
    }
    area /= 2.0;

    return area;

}

Point GeomCell::center() const
{
    Point ret;

    const size_t npoints = _boundary.size();
    for (size_t ipoint=0; ipoint < npoints; ++ipoint) {
	const Point& point = _boundary[ipoint];
	ret.x += point.x;
	ret.y += point.y;
	ret.z += point.z;
    }
    for (size_t ind=0; ind<3; ++ind) {
	ret.x /= npoints;
	ret.y /= npoints;
	ret.z /= npoints;
    }

    return ret;
}
