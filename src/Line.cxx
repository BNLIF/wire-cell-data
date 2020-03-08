#include "WCPData/Line.h"

using namespace WCP;

Line::Line(Point& p1, Point& p2)
  : p1(p1)
  , p2(p2)
{
  update_dir();
}

Line::Line(Point& p1, TVector3& dir)
  : p1(p1)
  , dir(dir)
{
  p2.x = p1.x +dir.X();
  p2.y = p1.y +dir.Y();
  p2.z = p1.z +dir.Z();
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

double Line::closest_dis(Line &l1){
  double dis = 0;
  
  TVector3 ca(p1.x-l1.get_p1().x,p1.y-l1.get_p1().y,p1.z-l1.get_p1().z);
  TVector3 bd = dir.Cross(l1.get_dir());
  dis = ca.Dot(bd)/bd.Mag();
  return fabs(dis);
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
  update_dir();
}

void Line::ReverseDir(){
  Point p3;
  p3 = p1;
  p1 = p2;
  p2 = p3;
  update_dir();
}
