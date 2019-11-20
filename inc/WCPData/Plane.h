#ifndef WCPData_Plane_h
#define WCPData_Plane_h

#include "WCPData/Line.h"

namespace WCP {
  class Plane{
  public:
    Plane(Point p1, Point p2, Point p3);
    
    bool sameline();
    
    Line& CrossLineCommonPoint(Plane& plane); 

    void cal_perp();
    TVector3& get_perp_vec(){return perp_vec;};

    Point& get_p1(){return p1;};
    Point& get_p2(){return p2;};
    Point& get_p3(){return p3;};

    ~Plane();
  protected:
    Point p1;
    Point p2;
    Point p3;
    
    TVector3 perp_vec;
  };
}

#endif
