#ifndef WIRECELL_PR3DCLUSTER_H
#define WIRECELL_PR3DCLUSTER_H 

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"

#include "TVector3.h"

namespace WireCell{
  class PR3DCluster{
  public:
    PR3DCluster(int cluster_id);
    ~PR3DCluster();
    void AddCell(SlimMergeGeomCell* mcell, int time_slice);
    // void AddCell(SlimMergeGeomCell* mcell, int *time_slices, int ntime_slice);
    int get_cluster_id(){return cluster_id;};
    int get_num_mcells(){return mcells.size();};
    int get_num_time_slices(){return time_cells_set_map.size();};

    void Remove_duplicated_mcells();
    SMGCSelection Is_Connected(PR3DCluster* cluster1, int offset);
    std::map<int,SMGCSet>& get_time_cells_set_map(){return time_cells_set_map;};
    SMGCSelection& get_mcells(){return mcells;};
    std::map<SlimMergeGeomCell*, std::set<int>>& get_cell_times_set_map(){return cell_times_set_map;};

    void Create_point_cloud();
    ToyPointCloud* get_point_cloud(){return point_cloud;};
    
    void Calc_PCA();
    Vector get_center(){return center;};
    Vector get_PCA_axis(int axis){return PCA_axis[axis];};

    //  std::pair<double,double> HoughTrans(Point& p, double dis);
    TVector3 calc_dir(Point& p_test, Point& p, double dis);
    std::pair<SlimMergeGeomCell*,Point> get_closest_point_mcell(Point& p_test);
    Point calc_ave_pos(Point& p, double dis);
    
  protected:
    int cluster_id;
    SMGCSelection mcells;
    /* std::map<int,SMGCSelection> time_cells_map; */
    /* std::map<SlimMergeGeomCell*, std::vector<int>> cell_times_map; */

    Vector center;
    Vector PCA_axis[3];
    
    ToyPointCloud *point_cloud;
    
    
    std::map<int,SMGCSet> time_cells_set_map;
    std::map<SlimMergeGeomCell*, std::set<int>> cell_times_set_map;   
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
