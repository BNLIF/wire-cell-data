#ifndef GeomWireCellData_Slim3DCluster_h
#define GeomWireCellData_Slim3DCluster_h

#include "WireCellData/Units.h"
#include "WireCellData/SlimMergeGeomCell.h"

#include <set>
#include <vector>

namespace WireCell {
  typedef std::vector<std::set<const SlimMergeGeomCell*>> SlimMergeCellCluster;

  class Slim3DCluster {
  public:
    Slim3DCluster(SlimMergeGeomCell& cell);
  
    ~Slim3DCluster();

    int AddCell(SlimMergeGeomCell &cell); // add a cell, 0, no need to add, 1 add in
    GeomCellSelection get_allcell();
    SlimMergeCellCluster get_ordercell(){return cluster;};

    void MergeCluster(Slim3DCluster& cluster1);
    

  protected:
    SlimMergeCellCluster cluster; // vector of time 
    GeomCellSelection gcluster; // all merged cell
  };

  typedef std::set<Slim3DCluster*> Slim3DClusterSet;
  typedef std::list<Slim3DCluster*> Slim3DClusterList;
};

#endif
