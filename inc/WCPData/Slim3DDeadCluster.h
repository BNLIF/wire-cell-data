#ifndef GeomWCPData_Slim3DDeadCluster_h
#define GeomWCPData_Slim3DDeadCluster_h

#include "WCPData/SlimMergeGeomCell.h"
#include <map>

namespace WCP{
  class Slim3DDeadCluster{
  public:
    Slim3DDeadCluster(SlimMergeGeomCell& cell, int time_slice);
    ~Slim3DDeadCluster();
    
    int AddCell(SlimMergeGeomCell &cell, int time_slice, int offset=1);
    void MergeCluster(Slim3DDeadCluster &cluster1);

    std::map<int,GeomCellSetp>& get_cluster(){return cluster;};
    GeomCellSetp& get_mcells(){return gcluster;};
    std::map<SlimMergeGeomCell*,std::set<int>>& get_mcell_time_map(){return mcell_time_map;};
    
    bool IsContain(SlimMergeGeomCell &cell, int time_slice);

    bool Extend(int time_slice);
    void set_id(int value){id = value;};
    int get_id(){return id;};
    
  protected:
    int id;
    GeomCellSetp gcluster;
    std::map<int,GeomCellSetp> cluster;
    std::map<SlimMergeGeomCell*,std::set<int>> mcell_time_map;
  };
  typedef std::set<Slim3DDeadCluster*> Slim3DDeadClusterSet;
  
}

#endif
