#ifndef GeomWireCellData_Projected2DCluster_h
#define GeomWireCellData_Projected2DCluster_h

#include "WireCellData/Units.h"
#include "WireCellData/SlimMergeGeomCell.h"

#include <vector>

namespace WireCell{
  class Projected2DCluster{
  public:
    Projected2DCluster(WirePlaneType_t plane_no, GeomCell *parent_cell);
    ~Projected2DCluster();
  
    void AddCell(SlimMergeGeomCell *mcell);
    
    GeomCell* GetParentCell(){return parent_cell;};
    WirePlaneType_t GetPlaneNo(){return plane_no;};
    
  protected:
    GeomCell *parent_cell;
    // which wire plane
    WirePlaneType_t plane_no;
    // range of time slices
    int time_slice_limit[2];
    // range of wires
    int wire_limit[2];

    std::vector<int> time_slice_array;
    std::vector<std::vector<std::pair<int,int>>> wire_array;
    
  };
}


#endif
