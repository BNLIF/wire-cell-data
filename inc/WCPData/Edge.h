#ifndef WCPData_Edge_h
#define WCPData_Edge_h

#include "WCPData/Point.h"

namespace WCP {
  typedef std::pair<WCP::Point,WCP::Point> Edge;
  typedef std::vector<WCP::Edge> EdgeVector;
  
  static bool operator==(const WCP::Edge& a, const WCP::Edge& b){
    if ((a.first==b.first && a.second == b.second) ||
    	(a.second == b.first && a.first == b.second)){
      return true;
    }else{
      return false;
    }
  }

  static bool CompareEdge(const WCP::Edge& a, const WCP::Edge& b){
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
  
  static bool ComparePoint(const WCP::Point& a, const WCP::Point& b){
    double dis = sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2)+pow(a.z-b.z,2))/units::m;
    if (dis < 0.0002){
      return true;
    }else{
      return false;
    }
  }

  
}

#endif
