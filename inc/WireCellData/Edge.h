#ifndef WireCellData_Edge_h
#define WireCellData_Edge_h

#include "WireCellData/Point.h"

namespace WireCell {
  typedef std::pair<WireCell::Point,WireCell::Point> Edge;
  typedef std::vector<WireCell::Edge> EdgeVector;
  
  static bool operator==(const WireCell::Edge& a, const WireCell::Edge& b){
    if ((a.first==b.first && a.second == b.second) ||
    	(a.second == b.first && a.first == b.second)){
      return true;
    }else{
      return false;
    }
  }
  static bool CompareEdge(const WireCell::Edge& a, const WireCell::Edge& b){
    double dis1,dis2,dis3,dis4;
    dis1 = sqrt(pow(a.first.x-b.first.x,2)+pow(a.first.y-b.first.y,2)+pow(a.first.z-b.first.z,2))/units::m;
    dis2 = sqrt(pow(a.second.x-b.second.x,2)+pow(a.second.y-b.second.y,2)+pow(a.second.z-b.second.z,2))/units::m;
    dis3 = sqrt(pow(a.first.x-b.second.x,2)+pow(a.first.y-b.second.y,2)+pow(a.first.z-b.second.z,2))/units::m;
    dis4 = sqrt(pow(a.second.x-b.first.x,2)+pow(a.second.y-b.first.y,2)+pow(a.second.z-b.first.z,2))/units::m;
    if ((dis1<0.0002 && dis2<0.0002) || (dis3 < 0.0002 && dis4 < 0.0002)){
      return true;
    }else{
      return false;
    }
  }
  
}

#endif
