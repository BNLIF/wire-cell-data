#ifndef GeomWireCellData_Slim3DDeadCluster_h
#define GeomWireCellData_Slim3DDeadCluster_h

#include "WireCellData/SlimMergeGeomCell.h"
#include <map>

namespace WireCell{
  class Slim3DDeadCluster{
  public:
    Slim3DDeadCluster(SlimMergeGeomCell& cell, int time_slice);
    ~Slim3DDeadCluster();
    
    int AddCell(SlimMergeGeomCell &cell, int time_slice, int offset=1);
    void MergeCluster(Slim3DDeadCluster &cluster1);
    
  protected:
    GeomCellSelection gcluster;
    std::map<int,GeomCellSelection> cluster;
  };
}

#endif
