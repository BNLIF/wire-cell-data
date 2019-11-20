#ifndef WCP_TOYPOINTCLOUD_H
#define WCP_TOYPOINTCLOUD_H


#include "WCPNanoflann/nanoflann.h" 
#include "WCPData/WCPointCloud.h"
#include "WCPQuickhull/QuickHull.h"
#include "WCPQuickhull/MathUtils.h"

#include "TVector3.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WCP::WCPointCloud<double> > ,
    WCP::WCPointCloud<double>,
    3 /* dim */
  > my_kd_tree_t;

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WCP::WC2DPointCloud<double> > ,
    WCP::WC2DPointCloud<double>,
    2 /* dim */
  > my_kd_tree_2d_t;

namespace WCP{
  
  class ToyPointCloud {
  public:
    ToyPointCloud(double angle_u = 1.0472, double angle_v = -1.0472, double angle_w = 0);
    ~ToyPointCloud();

    
    // build point cloud
    void AddPoint(WCPointCloud<double>::WCPoint& wcp, WC2DPointCloud<double>::WC2DPoint& wcp_u, WC2DPointCloud<double>::WC2DPoint& wcp_v, WC2DPointCloud<double>::WC2DPoint& wcp_w);
    void AddPoint(WCP::Point& p, std::tuple<int,int,int>& wires_index, WCP::SlimMergeGeomCell *mcell);
    void AddPoints(WCP::PointVector& ps, std::vector<std::tuple<int,int,int>>& wires_indices, WCP::SlimMergeGeomCell *mcell);
    void AddPoints(WCP::PointVector& ps);
    void build_kdtree_index();

    // find 3D oiubts
    std::tuple<int,int,double> get_closest_points(ToyPointCloud *point_could);
    std::pair<int,double> get_closest_point_along_vec(Point& p_test, TVector3& dir, double test_dis, double dis_step, double angle_cut, double dis_cut);
    
    
    // hull, not useful ...
    std::vector<WCPointCloud<double>::WCPoint> get_hull();

    // function to find 3D close points
    std::map<WCP::SlimMergeGeomCell*, WCP::Point> get_closest_mcell(WCP::Point& p, int N);
    std::map<WCP::SlimMergeGeomCell*, WCP::Point> get_closest_mcell(WCP::Point& p, double radius);
    std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> get_closest_points(WCP::Point& p, int N);
    std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> get_closest_points(WCP::Point& p, double radius);

    double get_closest_dis(WCP::Point& p);
    double get_closest_2d_dis(double x, double y, int plane);
    
    std::pair<double, WCP::Point> get_closest_point(WCP::Point& p);
    
    WCP::WCPointCloud<double>::WCPoint& get_closest_wcpoint(WCP::WCPointCloud<double>::WCPoint& wcp);
    WCP::WCPointCloud<double>::WCPoint& get_closest_wcpoint(WCP::Point& p);

    //function to find 2D points ... 
    std::vector<size_t> get_closest_2d_index(WCP::Point& p, int N, int plane);
    std::vector<size_t> get_closest_2d_index(WCP::Point& p, double radius, int plane);

    std::pair<int,double> get_closest_2d_dis(WCP::Point& p, int plane);
    
    int get_num_points(){return cloud.pts.size();};
    WCP::WCPointCloud<double>& get_cloud(){return cloud;};
    WCP::WC2DPointCloud<double>& get_cloud_u(){return cloud_u;};
    WCP::WC2DPointCloud<double>& get_cloud_v(){return cloud_v;};
    WCP::WC2DPointCloud<double>& get_cloud_w(){return cloud_w;};
    
    std::vector<int>& get_mcell_indices(WCP::SlimMergeGeomCell* mcell){
      return map_mcell_indices[mcell];
    };
    

    void Print();
    
  protected:
    std::vector<std::pair<size_t,double>> get_closest_index(WCP::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest_index(WCP::Point& p, double radius);

    
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, int N, int plane);
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, double radius, int plane);

    double angle_u, angle_v, angle_w;
    
    
    WCP::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
    WCP::WC2DPointCloud<double> cloud_u;
    my_kd_tree_2d_t *index_u;
    WCP::WC2DPointCloud<double> cloud_v;
    my_kd_tree_2d_t *index_v;
    WCP::WC2DPointCloud<double> cloud_w;
    my_kd_tree_2d_t *index_w;
    
    // map of mcells to WCPoints
    std::map<WCP::SlimMergeGeomCell*, std::vector<int>> map_mcell_indices; 
    
    

  };

  typedef std::vector<ToyPointCloud*> ToyPointCloudSelection;
   
}

#endif
