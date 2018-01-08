#ifndef WireCell_DYNAMICTOYPOINTCLOUD_H
#define WireCell_DYNAMICTOYPOINTCLOUD_H

#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellData/PR3DCluster.h"

#include "TVector3.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WCPointCloud<double> > ,
    WireCell::WCPointCloud<double>,
    3 /* dim */
  > my_dynamic_kd_tree_t;

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WC2DPointCloud<double> > ,
    WireCell::WC2DPointCloud<double>,
    2 /* dim */
  > my_dynamic_kd_tree_2d_t;


namespace WireCell{
  class DynamicToyPointCloud{
  public:
    DynamicToyPointCloud(double angle_u = 1.0472, double angle_v = -1.0472, double angle_w = 0);
    ~DynamicToyPointCloud();

  protected:
    double angle_u, angle_v, angle_w;
    WireCell::WCPointCloud<double> cloud;
    my_dynamic_kd_tree_t *index;
    WireCell::WC2DPointCloud<double> cloud_u;
    my_dynamic_kd_tree_2d_t *index_u;
    WireCell::WC2DPointCloud<double> cloud_v;
    my_dynamic_kd_tree_2d_t *index_v;
    WireCell::WC2DPointCloud<double> cloud_w;
    my_dynamic_kd_tree_2d_t *index_w;

    std::vector<PR3DCluster*> vec_index_cluster;
    
  };
}

#endif
