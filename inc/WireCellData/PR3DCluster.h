#ifndef WIRECELL_PR3DCLUSTER_H
#define WIRECELL_PR3DCLUSTER_H 

#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/ToyPointCloud.h"
#include "WireCellData/ToyCTPointCloud.h"
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
    void set_cluster_id(int value){cluster_id=value;};
    int get_num_mcells(){return mcells.size();};
    int get_num_points(){return point_cloud->get_num_points(); };
    int get_num_time_slices(){return time_cells_set_map.size();};

    bool judge_vertex(Point& p_test, double asy_cut = 1./3., double occupied_cut = 0.85);
    
    void Remove_duplicated_mcells();
    SMGCSelection Is_Connected(PR3DCluster* cluster1, int offset);
    std::map<int,SMGCSet>& get_time_cells_set_map(){return time_cells_set_map;};
    SMGCSelection& get_mcells(){return mcells;};
    std::map<SlimMergeGeomCell*, std::set<int>>& get_cell_times_set_map(){return cell_times_set_map;};

    void Create_point_cloud(WireCell::ToyPointCloud *global_point_cloud = 0);
    ToyPointCloud* get_point_cloud(){return point_cloud;};
    void Update_mcell_cluster_map(std::map<WireCell::SlimMergeGeomCell*,WireCell::PR3DCluster*>& mcell_cluster_map);
    
    void Create_graph();
    void Del_graph();
    void dijkstra_shortest_paths(WCPointCloud<double>::WCPoint& wcp_source);
    void cal_shortest_path(WCPointCloud<double>::WCPoint& wcp_target);
    void fine_tracking(int num_pts_cut = 10);
    
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> get_highest_lowest_wcps();
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> get_front_back_wcps();
    std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> get_earliest_latest_wcps();
    std::pair<Point,Point> get_two_extreme_points();
    
    
    
    WCPointCloud<double>::WCPoint get_furthest_wcpoint(WCPointCloud<double>::WCPoint wcp, TVector3 dir, double step = 5*units::cm, int allowed_nstep = 12);
    void adjust_wcpoints_parallel(WCPointCloud<double>::WCPoint& start_wcp, WCPointCloud<double>::WCPoint& end_wcp);

    bool Construct_skeleton();
    
    void Calc_PCA();
    bool flag_PCA;
    void Calc_PCA(WireCell::PointVector& points);
    Vector get_center(){return center;};
    Vector get_PCA_axis(int axis){return PCA_axis[axis];};
    double get_PCA_value(int axis){return PCA_values[axis];};

    std::pair<int,int> get_num_points(Point& p, TVector3& dir);
    std::pair<int,int> get_num_points(Point& p, TVector3& dir, double dis);
    
    std::pair<double,double> HoughTrans(Point& p, double dis);
    TVector3 VHoughTrans(Point& p, double dis);
    TVector3 calc_dir(Point& p_test, Point& p, double dis);
    TVector3 calc_PCA_dir(Point& p, double dis);

    std::pair<Point, double> get_closest_point_along_vec(Point& p_test, TVector3 dir, double test_dis, double dis_step, double angle_cut, double dis_cut);
    
    std::pair<SlimMergeGeomCell*,Point> get_closest_point_mcell(Point& p_test);
   
    
      
    Point calc_ave_pos(Point& p, double dis);
    Point calc_ave_pos(Point& p, int N);
    int get_num_points(Point& p_test, double dis);

    
    std::list<WCPointCloud<double>::WCPoint>& get_path_wcps(){return path_wcps;};
    std::list<SlimMergeGeomCell*>& get_path_mcells(){return path_mcells;};

    bool get_fine_tracking_flag(){return flag_fine_tracking;};
    PointVector& get_fine_tracking_path(){return fine_tracking_path;};
    TVector3 get_ft_dir_end(float mean_dis, float dis_cut);

    std::vector<PR3DCluster*> examine_x_boundary(double low_limit = -1*units::cm, double high_limit = 257*units::cm);
    
    std::vector<int> get_uvwt_range();

    void get_projection(std::vector<int>& proj_channel, std::vector<int>& proj_timeslice, std::vector<int>& proj_charge, std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map);
    
  protected:

    bool compare_edge(std::tuple<int,int,double>& index1, std::tuple<int,int,double>& index2);
    
    int cluster_id;
    
    SMGCSelection mcells;
    std::map<int,SMGCSet> time_cells_set_map;
    std::map<SlimMergeGeomCell*, std::set<int>> cell_times_set_map;   

    
    ToyPointCloud *point_cloud;
    Vector center;
    
    Vector PCA_axis[3];
    double PCA_values[3];
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
