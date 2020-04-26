#ifndef WCPData_Line_h
#define WCPData_Line_h

#include "WCPData/Point.h"
#include "TVector3.h"

namespace WCP {
  class Line {
  public:
    Line(Point& p1, Point& p2);
    Line(Point& p1, TVector3& dir);
    ~Line();

    void UpdatePoint(int flag, Point &p3);

    double closest_dis(Point &p3);
    double closest_dis(Line &l1);
    std::pair<Point, Point> closest_dis_points(Line &l1);
    
    TVector3& get_dir(){return dir;};
    Point& get_p1(){return p1;};
    Point& get_p2(){return p2;};
    
    void update_dir();
    TVector3& vec(){return dir;};
    void ReverseDir();

  protected:
    Point p1;
    Point p2;
    TVector3 dir;
  };
}


#endif
