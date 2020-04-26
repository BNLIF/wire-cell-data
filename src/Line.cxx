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

std::pair<Point, Point> Line::closest_dis_points(Line &l1){
  TVector3 d(p1.x-l1.get_p1().x,p1.y-l1.get_p1().y,p1.z-l1.get_p1().z);
  TVector3 c = dir.Cross(l1.get_dir()).Unit();
  TVector3 proj = d.Dot(l1.get_dir())/l1.get_dir().Mag() * l1.get_dir().Unit();
  TVector3 rej = d - proj - d.Dot(c)/c.Mag()*c.Unit();
  
  Point tp1;
  tp1.x = p1.x - rej.Mag() * dir.Unit().X()/dir.Unit().Dot(rej.Unit());
  tp1.y = p1.y - rej.Mag() * dir.Unit().Y()/dir.Unit().Dot(rej.Unit());
  tp1.z = p1.z - rej.Mag() * dir.Unit().Z()/dir.Unit().Dot(rej.Unit());

  TVector3 d_p = -d;
  TVector3 c_p = -c;
  TVector3 proj_p = d_p.Dot(dir)/dir.Mag() * dir.Unit();
  TVector3 rej_p = d_p - proj_p - d_p.Dot(c_p)/c_p.Mag()*c_p.Unit();
  
  Point tp2;
  tp2.x = l1.get_p1().x - rej_p.Mag() * l1.get_dir().Unit().X()/l1.get_dir().Unit().Dot(rej_p.Unit());
  tp2.y = l1.get_p1().y - rej_p.Mag() * l1.get_dir().Unit().Y()/l1.get_dir().Unit().Dot(rej_p.Unit());
  tp2.z = l1.get_p1().z - rej_p.Mag() * l1.get_dir().Unit().Z()/l1.get_dir().Unit().Dot(rej_p.Unit());

  // std::cout << sqrt(pow(tp1.x-tp2.x,2) + pow(tp1.y-tp2.y,2) + pow(tp1.z-tp2.z,2)) << std::endl;
  return std::make_pair(tp1, tp2);
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
