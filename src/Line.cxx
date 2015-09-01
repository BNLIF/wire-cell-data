#include "WireCellData/Line.h"

using namespace WireCell;

Line::Line(Point p1, Point p2)
  : p1(p1)
  , p2(p2)
{
  update_dir();
}

Line::~Line(){
};

void Line::update_dir(){
  dir.SetXYZ(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

double Line::closest_dis(Point &p3){
  double dis = 0;
  
  TVector3 d31(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
  TVector3 d32(p3.x-p2.x,p3.y-p2.y,p3.z-p2.z);

  TVector3 v = d31.Cross(d32);
  
  dis = v.Mag()/dir.Mag();
  
  return dis;
}


void Line::UpdatePoint(int flag, Point &p3){
  if (flag == 1){
    p1.x = p3.x;
    p1.y = p3.y;
    p1.z = p3.z;
  }else if (flag == 2){
    p2.x = p3.x;
    p2.y = p3.y;
    p2.z = p3.z;
  }
}

void Line::ReverseDir(){
  Point p3;
  p3 = p1;
  p1 = p2;
  p2 = p3;
}
