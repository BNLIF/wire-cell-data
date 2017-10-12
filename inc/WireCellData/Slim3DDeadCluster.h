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

    std::map<int,GeomCellSetp> get_cluster(){return cluster;};
    GeomCellSetp get_mcells(){return gcluster;};
    std::map<SlimMergeGeomCell*,std::set<int>> get_mcell_time_map(){return mcell_time_map;};
    
  protected:
    GeomCellSetp gcluster;
    std::map<int,GeomCellSetp> cluster;
    std::map<SlimMergeGeomCell*,std::set<int>> mcell_time_map;
  };
  typedef std::set<Slim3DDeadCluster*> Slim3DDeadClusterSet;
  
}

#endif
