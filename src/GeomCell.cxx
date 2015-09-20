#include "WireCellData/GeomCell.h"

#include <vector>
#include <cmath>

using namespace std;
using namespace WireCell;

GeomCell::GeomCell(int ident, const PointVector& boundary, int flag)
    : _ident(ident)
    , _boundary(boundary)
{
  order_boundary();
  if (flag == 1){
    if (boundary.size()>2){
      for (int i=0;i<boundary.size()-1;i++){
	Edge a(_boundary.at(i),_boundary.at(i+1));
	_edge.push_back(a);
      }
      Edge a(_boundary.at(boundary.size()-1),_boundary.at(0));
      _edge.push_back(a);
    }
  }
  flag_center = 0;
  flag_cross_section = 0;
  ret.x = 0;
  ret.y = 0;
  ret.z = 0;
  area = 0;

  uwire = 0;
  vwire = 0;
  wwire = 0;
}


GeomCell::GeomCell(const GeomCell *cell, int flag){
  _ident = cell->ident();
  _boundary = cell->boundary();
  order_boundary();
  if (flag == 1){
    if (_boundary.size()>2){
      for (int i=0;i<_boundary.size()-1;i++){
	Edge a(_boundary.at(i),_boundary.at(i+1));
	_edge.push_back(a);
      }
      Edge a(_boundary.at(_boundary.size()-1),_boundary.at(0));
      _edge.push_back(a);
    }
  }
  flag_center = 0;
  flag_cross_section = 0;
  ret.x = 0;
  ret.y = 0;
  ret.z = 0;
  area = 0;

  uwire = 0;
  vwire = 0;
  wwire = 0;
}

GeomCell::~GeomCell()
{
  _boundary.clear();
  _edge.clear();
}

std::ostream & WireCell::operator<<(std::ostream &os, const GeomCell& gc)
{
    return os << "<WireCell::GeomCell " << gc.ident() << ">";
}

double GeomCell::cross_section() const
{
  // double area = 0.0;
  if (flag_cross_section == 0){
    const size_t npoints = _boundary.size();
    int prev = npoints - 1;
    
    for (int ind = 0; ind < npoints; ++ind) {
      double z = _boundary.at(prev).z + _boundary.at(ind).z;
      double y = _boundary.at(prev).y - _boundary.at(ind).y;
      area += y*z;
      prev = ind;
    }
    area /= 2.0;
    flag_cross_section = 1;
  }
  return fabs(area);

}

Point GeomCell::center() const
{
  //Point ret(0,0,0);
  if (flag_center ==0){
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
    flag_center = 1;
  }
  
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
