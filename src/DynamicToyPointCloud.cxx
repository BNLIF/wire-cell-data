#include "WireCellData/DynamicToyPointCloud.h"

using namespace WireCell;

WireCell::DynamicToyPointCloud::DynamicToyPointCloud(double angle_u, double angle_v, double angle_w)
  : angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
{
  index = new my_dynamic_kd_tree_t(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_u = new my_dynamic_kd_tree_2d_t(2, cloud_u,  nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_v = new my_dynamic_kd_tree_2d_t(2, cloud_v,  nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_w = new my_dynamic_kd_tree_2d_t(2, cloud_w,  nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
      
}

WireCell::DynamicToyPointCloud::~DynamicToyPointCloud(){
  delete index;
  delete index_u;
  delete index_v;
  delete index_w;
  
  cloud.pts.clear();
  cloud_u.pts.clear();
  cloud_v.pts.clear();
  cloud_w.pts.clear();
}
