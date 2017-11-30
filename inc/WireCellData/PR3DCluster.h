#ifndef WIRECELL_PR3DCLUSTER_H
#define WIRECELL_PR3DCLUSTER_H 

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"
#include "TVector3.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;



namespace WireCell{

  struct VertexProp {
    int index;
    //WCPointCloud<double>::WCPoint wcpoint;
    // add pointer to merged cell
  };
  struct EdgeProp {
    float dist; // edge distance
  };
  
  typedef adjacency_list<vecS, vecS, undirectedS, VertexProp, EdgeProp> MCUGraph;
  typedef graph_traits<MCUGraph>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<MCUGraph>::edge_descriptor edge_descriptor;
  
  
  
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

    void Create_graph();
    void dijkstra_shortest_paths(WCPointCloud<double>::WCPoint& wcp_source);
    void cal_shortest_path(WCPointCloud<double>::WCPoint& wcp_target);
    void fine_tracking();
    
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> get_highest_lowest_wcps();
    
    void Calc_PCA();
    Vector get_center(){return center;};
    Vector get_PCA_axis(int axis){return PCA_axis[axis];};

    std::pair<double,double> HoughTrans(Point& p, double dis);
    TVector3 VHoughTrans(Point& p, double dis);
    TVector3 calc_dir(Point& p_test, Point& p, double dis);
    TVector3 calc_PCA_dir(Point& p, double dis);
    
    std::pair<SlimMergeGeomCell*,Point> get_closest_point_mcell(Point& p_test);
    Point calc_ave_pos(Point& p, double dis);

    std::list<WCPointCloud<double>::WCPoint>& get_path_wcps(){return path_wcps;};
    std::list<SlimMergeGeomCell*>& get_path_mcells(){return path_mcells;};

    bool get_fine_tracking_flag(){return flag_fine_tracking;};
    PointVector& get_fine_tracking_path(){return fine_tracking_path;};

    std::vector<int> get_uvwt_range();
    
  protected:
    
    int cluster_id;
    
    SMGCSelection mcells;
    std::map<int,SMGCSet> time_cells_set_map;
    std::map<SlimMergeGeomCell*, std::set<int>> cell_times_set_map;   

    
    ToyPointCloud *point_cloud;
    Vector center;
    
    Vector PCA_axis[3];

    // graph 
    MCUGraph *graph;

    // create things for Dijkstra
    std::vector<vertex_descriptor> parents;
    std::vector<int> distances;
    int source_wcp_index;
    // return ... 
    std::list<WCPointCloud<double>::WCPoint> path_wcps;
    std::list<SlimMergeGeomCell*> path_mcells;
    int dest_wcp_index;
    // fine tracking related ...
    bool flag_fine_tracking;
    PointVector fine_tracking_path;

    
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}

#endif
