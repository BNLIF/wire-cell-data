#ifndef GeomWireCellData_Projected2DCluster_h
#define GeomWireCellData_Projected2DCluster_h

#include "WireCellData/Units.h"
#include "WireCellData/SlimMergeGeomCell.h"

#include <vector>
#include <list>
#include <map>

namespace WireCell{
  class Projected2DCluster{
  public:
    Projected2DCluster(WirePlaneType_t plane_no);
    ~Projected2DCluster();
  
    void AddCell(SlimMergeGeomCell *mcell);

    // +1 the 
    int judge_coverage(Projected2DCluster *cluster);
    
    
    WirePlaneType_t GetPlaneNo(){return plane_no;};

    int get_number_time_slices(){return time_slice_array.size();};
    int get_low_time_limit(){return time_slice_limit[0];}
    int get_high_time_limit(){return time_slice_limit[1];}
    int get_low_wire_limit(){return wire_limit[0];};
    int get_high_wire_limit(){return wire_limit[1];};
    std::vector<int>& get_time_slice_array(){return time_slice_array;};
    std::map<int, std::list<std::pair<int,int>>>& get_time_wires_map(){return time_wires_map;};
    
    void Print();
    
    
  protected:
    GeomCell *parent_cell;
    // which wire plane
    WirePlaneType_t plane_no;
    // range of time slices
    int time_slice_limit[2];
    // range of wires
    int wire_limit[2];

    std::vector<int> time_slice_array;
    std::map<int,std::list<std::pair<int,int>>> time_wires_map;
    
  };
}


#endif
