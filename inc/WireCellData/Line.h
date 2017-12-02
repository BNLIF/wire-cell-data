#ifndef WireCellData_Line_h
#define WireCellData_Line_h

#include "WireCellData/Point.h"
#include "TVector3.h"

namespace WireCell {
  class Line {
  public:
    Line(Point& p1, Point& p2);
    Line(Point& p1, TVector3& dir);
    ~Line();

    void UpdatePoint(int flag, Point &p3);

    double closest_dis(Point &p3);
    double closest_dis(Line &l1);
    
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
