#ifndef GeomWCPData_Slim3DDeadCluster_h
#define GeomWCPData_Slim3DDeadCluster_h

#include "WCPData/SlimMergeGeomCell.h"
#include <map>

namespace WCP{
  class Slim3DDeadCluster{
  public:
    Slim3DDeadCluster(int ident, SlimMergeGeomCell& cell, int time_slice);
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
    int get_ident() const{return ident;};
    
  protected:
    int id;
    int ident;
    GeomCellSetp gcluster;
    std::map<int,GeomCellSetp> cluster;
    std::map<SlimMergeGeomCell*,std::set<int>> mcell_time_map;
  };

   // ensure the order ...
    struct Slim3DDeadClusterComparep {
      bool operator() (const Slim3DDeadCluster* a, const Slim3DDeadCluster* b) const {

	if (a && b){
	  if (a->get_ident() != b->get_ident())
	    return a->get_ident() < b->get_ident();
	  else
	    return a < b;
	}
	return false;
      }
    };
  
    typedef std::set<Slim3DDeadCluster*, Slim3DDeadClusterComparep> Slim3DDeadClusterSet;
  
}

#endif
