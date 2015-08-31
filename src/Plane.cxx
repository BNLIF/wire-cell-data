#include "WireCellData/Plane.h"

using namespace WireCell;

Plane::Plane(Point &p1, Point &p2, Point &p3)
  : p1(p1)
  , p2(p2)
  , p3(p3)
{
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
