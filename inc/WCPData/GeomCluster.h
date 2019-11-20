#ifndef GeomWCPData_Cluster_h
#define GeomWCPData_Cluster_h

#include "WCPData/Units.h"
#include "WCPData/MergeGeomCell.h"

#include <set>
#include <vector>

namespace WCP {
  typedef std::vector<std::set<const MergeGeomCell*>> MergeCellCluster;

  class GeomCluster {
  public:
    GeomCluster(MergeGeomCell& cell);
  
    ~GeomCluster();

    int AddCell(MergeGeomCell &cell); // add a cell, 0, no need to add, 1 add in
    GeomCellSelection get_allcell();
    MergeCellCluster get_ordercell(){return cluster;};

    void MergeCluster(GeomCluster& cluster1);
    

  protected:
    MergeCellCluster cluster; // vector of time 
    GeomCellSelection gcluster; // all merged cell
  };

  typedef std::set<GeomCluster*> GeomClusterSet;
  typedef std::list<GeomCluster*> GeomClusterList;
};

#endif
