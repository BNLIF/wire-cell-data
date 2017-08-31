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
    Projected2DCluster(WirePlaneType_t plane_no, int id);
    ~Projected2DCluster();
  
    void AddCell(SlimMergeGeomCell *mcell);

    // +1 the 
    int judge_coverage(Projected2DCluster *cluster);

    std::vector<int> calc_coverage(Projected2DCluster *cluster);
    int judge_coverage_alt(Projected2DCluster *cluster);
    
    WirePlaneType_t GetPlaneNo(){return plane_no;};

    int get_number_time_slices(){return time_slice_array.size();};
    int get_low_time_limit(){return time_slice_limit[0];}
    int get_high_time_limit(){return time_slice_limit[1];}
    int get_low_wire_limit(){return wire_limit[0];};
    int get_high_wire_limit(){return wire_limit[1];};
    std::vector<int>& get_time_slice_array(){return time_slice_array;};
    std::map<int, std::list<std::pair<int,int>>>& get_time_wires_map(){return time_wires_map;};
    int get_charge(std::pair<int,int> time_wire){return time_wire_charge_map[time_wire];};
    void Print();
    int get_parent_cluster_id(){return parent_cluster_id;};
    void set_estimated_total_charge(float value){estimated_total_charge = value;};
    float get_estimated_total_charge(){return estimated_total_charge;};
    int calc_dead_time_wire();
    
  protected:
    float estimated_total_charge;
    int parent_cluster_id;
    //GeomCell *parent_cell;
    // which wire plane
    WirePlaneType_t plane_no;
    // range of time slices
    int time_slice_limit[2];
    // range of wires
    int wire_limit[2];

    std::vector<int> time_slice_array;
    std::map<int,std::list<std::pair<int,int>>> time_wires_map;

    std::vector<int> dead_time_slice_array;
    std::map<int,std::list<std::pair<int,int>>> dead_time_wires_map;
    
    std::map<std::pair<int,int>,float> time_wire_charge_map;
  };
}


#endif
