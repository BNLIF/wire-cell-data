#ifndef WireCell_TOYCTPOINTCLOUD_H
#define WireCell_TOYCTPOINTCLOUD_H

#include "WireCellNanoflann/nanoflann.h"
#include "WireCellData/WCPointCloud.h"

#include <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
nanoflann::L2_Simple_Adaptor<double, WireCell::CTPointCloud<double> > ,
  WireCell::CTPointCloud<double>,
  2 /* dim */
  > my_kd_tree_ct_t;

namespace WireCell{
  class ToyCTPointCloud {
  public:
    ToyCTPointCloud(int u_min_ch, int u_max_ch, int v_min_ch, int v_max_ch, int w_min_ch, int w_max_ch, double pitch_u, double pitch_v, double pitch_w, double offset_t, double convert_t);
    ~ToyCTPointCloud();

  protected:
    int u_min_ch, u_max_ch;
    int v_min_ch, v_max_ch;
    int w_min_ch, w_max_ch;
    double pitch_u, pitch_v, pitch_w;
    double offset_t, convert_t;

    
    WireCell::CTPointCloud<double> cloud_u;
    my_kd_tree_ct_t *index_u;
    WireCell::CTPointCloud<double> cloud_v;
    my_kd_tree_ct_t *index_v;
    WireCell::CTPointCloud<double> cloud_w;
    my_kd_tree_ct_t *index_w;
  };
}


#endif

