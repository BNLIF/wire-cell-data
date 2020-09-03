#ifndef GeomWCPData_Slim3DCluster_h
#define GeomWCPData_Slim3DCluster_h

#include "WCPData/Units.h"
#include "WCPData/SlimMergeGeomCell.h"
#include "WCPData/Projected2DCluster.h"
#include "WCPData/Slim3DDeadCluster.h"

#include <set>
#include <vector>

namespace WCP {
  typedef std::vector<SMGCSet> SlimMergeCellCluster;

  class Slim3DCluster {
  public:
    Slim3DCluster(int ident, SlimMergeGeomCell& cell);
  
    ~Slim3DCluster();

    int AddCell(SlimMergeGeomCell &cell, int offset=1); // add a cell, 0, no need to add, 1 add in
    void DirectAddCell(SlimMergeGeomCell &cell);
    void DirectOrderCells();
    
    GeomCellSelection& get_allcell();
    SlimMergeCellCluster& get_ordercell(){return cluster;};

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
    int get_ident() const{return ident;};
    
    float get_total_charge(){return total_charge;};
    float get_min_total_charge(){return min_total_charge;};

    void Form_maps(int offset, std::map<const GeomCell*, GeomCellSelection, WCP::GeomCellComparep>& front_cell_map, std::map<const GeomCell*, GeomCellSelection, WCP::GeomCellComparep>& back_cell_map);

    GeomCellSelection Is_Connected(Slim3DDeadCluster* cluster1 , int offset=1);
    
  protected:
    int id;

    int ident;
    
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

    /* std::map<const GeomCell*, GeomCellSelection> front_cell_map; */
    /* std::map<const GeomCell*, GeomCellSelection> back_cell_map; */
    
  };

   struct Slim3DClusterComparep {
      bool operator() (const Slim3DCluster* a, const Slim3DCluster* b) const {

	if (a && b){
	  if (a->get_ident() != b->get_ident())
	    return a->get_ident() < b->get_ident();
	  else
	    return a < b;
	}
	return false;
      }
    };

  
   typedef std::set<Slim3DCluster*, Slim3DClusterComparep> Slim3DClusterSet;
   typedef std::list<Slim3DCluster*> Slim3DClusterList;
};

#endif
