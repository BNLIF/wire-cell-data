#include "WireCellData/ToyCTPointCloud.h"

using namespace WireCell;

ToyCTPointCloud::ToyCTPointCloud(int u_min_ch, int u_max_ch, int v_min_ch, int v_max_ch, int w_min_ch, int w_max_ch, double pitch_u, double pitch_v, double pitch_w, double offset_t, double convert_t)
  : u_min_ch(u_min_ch)
  , u_max_ch(u_max_ch)
  , v_min_ch(v_min_ch)
  , v_max_ch(v_max_ch)
  , w_min_ch(w_min_ch)
  , w_max_ch(w_max_ch)
  , pitch_u(pitch_u)
  , pitch_v(pitch_v)
  , pitch_w(pitch_w)
  , offset_t(offset_t)
  , convert_t(convert_t)
{
  index_u = 0;
  index_v = 0;
  index_w = 0;
}

ToyCTPointCloud::~ToyCTPointCloud(){
  if (index_u!=0)    delete index_u;
  if (index_v!=0)    delete index_v;
  if (index_w!=0)    delete index_w;

  cloud_u.pts.clear();
  cloud_v.pts.clear();
  cloud_w.pts.clear();
}
