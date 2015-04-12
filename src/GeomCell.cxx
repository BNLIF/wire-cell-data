#include "WireCellData/GeomCell.h"

#include <vector>
#include <cmath>

using namespace std;
using namespace WireCell;

GeomCell::GeomCell(int ident, const PointVector& boundary)
    : _ident(ident)
    , _boundary(boundary)
{
  order_boundary();
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

    return fabs(area);

}

Point GeomCell::center() const
{
  Point ret(0,0,0);
  
  const size_t npoints = _boundary.size();
  for (size_t ipoint=0; ipoint < npoints; ++ipoint) {
    const Point& point = _boundary[ipoint];
    ret.x += point.x;
    ret.y += point.y;
    ret.z += point.z;
    //std::cout << "qx1 " << point.y << " " << ret.y << std::endl;
  }
 
  ret.x /= npoints;
  ret.y /= npoints;
  ret.z /= npoints;
  
  
  return ret;
}

int GeomCell::order_boundary(){
  Point Center = center();
  
  std::map<float,Point> phi_boundary;

  for (int i=0;i!=_boundary.size();i++){
    Point p = _boundary[i];
    float phi = std::atan2(p.z-Center.z,p.y-Center.y);
    phi_boundary[phi] = p;
  }

  int i=0;
  for (std::map<float,Point>::iterator it = phi_boundary.begin(); it != phi_boundary.end(); ++it){
    _boundary[i] = it->second;
    i++;
  }

  return 0;
}


//  return std::atan2(y, x);
