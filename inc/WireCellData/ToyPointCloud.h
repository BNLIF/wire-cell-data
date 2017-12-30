#ifndef WireCell_TOYPOINTCLOUD_H
#define WireCell_TOYPOINTCLOUD_H


#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WCPointCloud<double> > ,
    WireCell::WCPointCloud<double>,
    3 /* dim */
  > my_kd_tree_t;

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WC2DPointCloud<double> > ,
    WireCell::WC2DPointCloud<double>,
    2 /* dim */
  > my_kd_tree_2d_t;

namespace WireCell{
  
  class ToyPointCloud {
  public:
    ToyPointCloud(double angle_u = 1.0472, double angle_v = -1.0472, double angle_w = 0);
    ~ToyPointCloud();

    
    // build point cloud
    void AddPoint(WCPointCloud<double>::WCPoint& wcp, WC2DPointCloud<double>::WC2DPoint& wcp_u, WC2DPointCloud<double>::WC2DPoint& wcp_v, WC2DPointCloud<double>::WC2DPoint& wcp_w);
    void AddPoint(WireCell::Point& p, std::tuple<int,int,int>& wires_index, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps, std::vector<std::tuple<int,int,int>>& wires_indices, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps);
    void build_kdtree_index();

    // find 3D oiubts
    std::tuple<int,int,double> get_closest_points(ToyPointCloud *point_could);

    // hull, not useful ...
    std::vector<WCPointCloud<double>::WCPoint> get_hull();

    // function to find 3D close points
    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, int N);
    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, double radius);
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, int N);
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, double radius);
    WireCell::WCPointCloud<double>::WCPoint& get_closest_wcpoint(WireCell::WCPointCloud<double>::WCPoint& wcp);
    WireCell::WCPointCloud<double>::WCPoint& get_closest_wcpoint(WireCell::Point& p);

    //function to find 2D points ... 
    std::vector<size_t> get_closest_2d_index(WireCell::Point& p, int N, int plane);
    std::vector<size_t> get_closest_2d_index(WireCell::Point& p, double radius, int plane);

    std::pair<int,double> get_closest_2d_dis(WireCell::Point& p, int plane);
    
    int get_num_points(){return cloud.pts.size();};
    WireCell::WCPointCloud<double>& get_cloud(){return cloud;};
    WireCell::WC2DPointCloud<double>& get_cloud_u(){return cloud_u;};
    WireCell::WC2DPointCloud<double>& get_cloud_v(){return cloud_v;};
    WireCell::WC2DPointCloud<double>& get_cloud_w(){return cloud_w;};
    
    std::vector<int>& get_mcell_indices(WireCell::SlimMergeGeomCell* mcell){
      return map_mcell_indices[mcell];
    };
    

    void Print();
    
  protected:
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, double radius);

    
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, int N, int plane);
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, double radius, int plane);

    double angle_u, angle_v, angle_w;
    
    
    WireCell::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
    WireCell::WC2DPointCloud<double> cloud_u;
    my_kd_tree_2d_t *index_u;
    WireCell::WC2DPointCloud<double> cloud_v;
    my_kd_tree_2d_t *index_v;
    WireCell::WC2DPointCloud<double> cloud_w;
    my_kd_tree_2d_t *index_w;
    
    // map of mcells to WCPoints
    std::map<WireCell::SlimMergeGeomCell*, std::vector<int>> map_mcell_indices; 
    
    

  };

  typedef std::vector<ToyPointCloud*> ToyPointCloudSelection;
   
}

#endif
