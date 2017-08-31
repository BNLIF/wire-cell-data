#ifndef GeomWireCellData_Slim3DCluster_h
#define GeomWireCellData_Slim3DCluster_h

#include "WireCellData/Units.h"
#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/Projected2DCluster.h"

#include <set>
#include <vector>

namespace WireCell {
  typedef std::vector<std::set<SlimMergeGeomCell*>> SlimMergeCellCluster;

  class Slim3DCluster {
  public:
    Slim3DCluster(SlimMergeGeomCell& cell);
  
    ~Slim3DCluster();

    int AddCell(SlimMergeGeomCell &cell, int offset=1); // add a cell, 0, no need to add, 1 add in
    GeomCellSelection get_allcell();
    SlimMergeCellCluster get_ordercell(){return cluster;};

    void MergeCluster(Slim3DCluster& cluster1);

    void Calc_Projection();
    Projected2DCluster* get_projection(WirePlaneType_t plane);
    int get_flag_saved(){return flag_saved;};
    void set_flag_saved(int value){flag_saved = value;};

    int get_flag_saved_1(){return flag_saved_1;};
    void set_flag_saved_1(int value){flag_saved_1 = value;};

    int get_flag_saved_u(){return flag_saved_u;};
    void set_flag_saved_u(int value){flag_saved_u = value;};
    int get_flag_saved_v(){return flag_saved_v;};
    void set_flag_saved_v(int value){flag_saved_v = value;};
    int get_flag_saved_w(){return flag_saved_w;};
    void set_flag_saved_w(int value){flag_saved_w = value;};
    void set_id(int value){id = value;};
    int get_id(){return id;};

    float get_total_charge(){return total_charge;};
    float get_min_total_charge(){return min_total_charge;};
    
  protected:
    int id;
    
    int flag_saved;
    int flag_saved_1;
    int flag_saved_u;
    int flag_saved_v;
    int flag_saved_w;


    float total_charge;
    float min_total_charge;
    
    SlimMergeCellCluster cluster; // vector of time 
    GeomCellSelection gcluster; // all merged cell

    Projected2DCluster *u_proj;
    Projected2DCluster *v_proj;
    Projected2DCluster *w_proj;
    
  };

  typedef std::set<Slim3DCluster*> Slim3DClusterSet;
  typedef std::list<Slim3DCluster*> Slim3DClusterList;
};

#endif
