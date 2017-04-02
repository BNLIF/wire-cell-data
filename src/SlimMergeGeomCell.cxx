#include "WireCellData/SlimMergeGeomCell.h"

using namespace WireCell;

WireCell::SlimMergeGeomCell::SlimMergeGeomCell(){
  
}

WireCell::SlimMergeGeomCell::~SlimMergeGeomCell(){
  uwires.clear();
  vwires.clear();
  wwires.clear();
}

void WireCell::SlimMergeGeomCell::AddBoundary(const PointVector& boundary){
  _boundary = boundary;
  flag_center = 0;
  flag_cross_section = 0;
  
}

// Point WireCell::SlimMergeGeomCell::center() const
// {
//   //Point ret(0,0,0);
//   if (flag_center ==0){
//     const size_t npoints = _boundary.size();
//     for (size_t ipoint=0; ipoint < npoints; ++ipoint) {
//       const Point& point = _boundary[ipoint];
//       ret.x += point.x;
//       ret.y += point.y;
//       ret.z += point.z;
//       //std::cout << "qx1 " << point.y << " " << ret.y << std::endl;
//     }
    
//     ret.x /= npoints;
//     ret.y /= npoints;
//     ret.z /= npoints;
//     flag_center = 1;
//   }
  
//   return ret;
// }

// int WireCell::SlimMergeGeomCell::order_boundary(){
//   Point Center = center();
  
//   std::map<float,Point> phi_boundary;

//   for (int i=0;i!=_boundary.size();i++){
//     Point p = _boundary[i];
//     float phi = std::atan2(p.z-Center.z,p.y-Center.y);
//     phi_boundary[phi] = p;
//   }

//   int i=0;
//   for (std::map<float,Point>::iterator it = phi_boundary.begin(); it != phi_boundary.end(); ++it){
//     _boundary[i] = it->second;
//     i++;
//   }

//   return 0;
// }


void WireCell::SlimMergeGeomCell::AddWire(const GeomWire *wire, WirePlaneType_t plane){
  if (plane == WirePlaneType_t(0)){
    if (find(uwires.begin(),uwires.end(),wire)==uwires.end())
      uwires.push_back(wire);
  }else if (plane == WirePlaneType_t(1)){
    if (find(vwires.begin(),vwires.end(),wire)==vwires.end())
      vwires.push_back(wire);
  }else if (plane == WirePlaneType_t(2)){
    if (find(wwires.begin(),wwires.end(),wire)==wwires.end())
      wwires.push_back(wire);
  }
}


void WireCell::SlimMergeGeomCell::OrderWires(){
  WireCell::sort_by_ident(uwires);
  WireCell::sort_by_ident(vwires);
  WireCell::sort_by_ident(wwires);
}
