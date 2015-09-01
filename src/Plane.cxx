#include "WireCellData/Plane.h"

using namespace WireCell;

Plane::Plane(Point p1, Point p2, Point p3)
  : p1(p1)
  , p2(p2)
  , p3(p3)
{
  cal_perp();
}

void Plane::cal_perp(){
  TVector3 v31(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
  TVector3 v32(p3.x-p2.x,p3.y-p2.y,p3.z-p2.z);
  
  perp_vec = v31.Cross(v32);
}

Plane::~Plane(){
}

bool Plane::sameline(){
  TVector3 v31(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
  TVector3 v32(p3.x-p2.x,p3.y-p2.y,p3.z-p2.z);

  TVector3 v = v31.Cross(v32);
  if (v.Mag()==0){
    return true;
  }else{
    return false;
  }
}


Line& Plane::CrossLineCommonPoint(Plane& plane1){
  // find the common points among two planes, 
  Point& pp1 = plane1.get_p1();
  Point& pp2 = plane1.get_p2();
  Point& pp3 = plane1.get_p3();

  Point common_point;
  
  if (pp1.x == p1.x && pp1.y == p1.y && pp1.z==p1.z){
    common_point = pp1;
  }else if (pp1.x == p2.x && pp1.y == p2.y && pp1.z==p2.z){
    common_point = pp1;
  }else if (pp1.x == p3.x && pp1.y == p3.y && pp1.z==p3.z){
    common_point = pp1;
  }else if (pp2.x == p1.x && pp2.y == p1.y && pp2.z==p1.z){
    common_point = pp2;
  }else if (pp2.x == p2.x && pp2.y == p2.y && pp2.z==p2.z){
    common_point = pp2;
  }else if (pp2.x == p3.x && pp2.y == p3.y && pp2.z==p3.z){
    common_point = pp2;
  }else if (pp3.x == p1.x && pp3.y == p1.y && pp3.z==p1.z){
    common_point = pp3;
  }else if (pp3.x == p2.x && pp3.y == p2.y && pp3.z==p2.z){
    common_point = pp3;
  }else if (pp3.x == p3.x && pp3.y == p3.y && pp3.z==p3.z){
    common_point = pp3;
  }

  // define the two perpendicular vector, and calcualte the perpendicular vector
  TVector3& perp_vec1 = plane1.get_perp_vec();
  TVector3 v = perp_vec.Cross(perp_vec1);

  Point new_point;
  new_point.x = v.x() + common_point.x;
  new_point.y = v.y() + common_point.y;
  new_point.z = v.z() + common_point.z;
  
  Line *line = new Line(common_point, new_point);
  return *line;
}
