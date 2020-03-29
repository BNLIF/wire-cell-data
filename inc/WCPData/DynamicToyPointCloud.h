#ifndef WCP_DYNAMICTOYPOINTCLOUD_H
#define WCP_DYNAMICTOYPOINTCLOUD_H

#include "WCPNanoflann/nanoflann.h" 
#include "WCPData/WCPointCloud.h"
#include "WCPData/PR3DCluster.h"

#include "TVector3.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WCP::WCPointCloud<double> > ,
    WCP::WCPointCloud<double>,
    3 /* dim */
  > my_dynamic_kd_tree_t;

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WCP::WC2DPointCloud<double> > ,
    WCP::WC2DPointCloud<double>,
    2 /* dim */
  > my_dynamic_kd_tree_2d_t;


namespace WCP{
  class DynamicToyPointCloud{
  public:
    DynamicToyPointCloud(double angle_u = 1.0472, double angle_v = -1.0472, double angle_w = 0);
    ~DynamicToyPointCloud();

    void AddPoints(PR3DCluster* cluster, int flag=0, double step = 0.6*units::cm); // flag 1 points, flag 2 scheleton
    void AddPoints(PointVector& pts, double step = 0.6*units::cm); // directly add points ...
    void AddPoints(ToyPointCloud* point_cloud, double step = 0.6*units::cm);
    
    void AddPoints(PR3DCluster* cluster, Point& p_test, TVector3& dir, double range, double step, double angle);

    std::tuple<double, PR3DCluster*, size_t> get_closest_point_info(WCP::Point& p);
    std::tuple<double, PR3DCluster*, size_t> get_closest_2d_point_info(WCP::Point& p, int plane);
    std::vector<std::tuple<double, PR3DCluster*, size_t>> get_2d_points_info(WCP::Point& p, double radius,  int plane);
    
    
    PR3DCluster* get_cluster(int index);

    std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> get_closest_points(WCP::Point& p, double radius);
    
    
    std::pair<double,double> HoughTrans(Point& p, double dis);
    TVector3 VHoughTrans(Point& p, double dis);

    WCP::WCPointCloud<double>& get_cloud(){return cloud;};
    
  protected:
    std::vector<std::pair<size_t,double>> get_closest_index(WCP::Point& p, int N); 
    std::vector<std::pair<size_t,double>> get_closest_index(WCP::Point& p, double radius); 

    
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, int N, int plane); 
    std::vector<std::pair<size_t,double>> get_closest_2d_index(double x, double y, double radius, int plane);

    
    double angle_u, angle_v, angle_w;
    WCP::WCPointCloud<double> cloud;
    my_dynamic_kd_tree_t *index;
    WCP::WC2DPointCloud<double> cloud_u;
    my_dynamic_kd_tree_2d_t *index_u;
    WCP::WC2DPointCloud<double> cloud_v;
    my_dynamic_kd_tree_2d_t *index_v;
    WCP::WC2DPointCloud<double> cloud_w;
    my_dynamic_kd_tree_2d_t *index_w;

    std::vector<PR3DCluster*> vec_index_cluster;
    
  };
}

#endif
