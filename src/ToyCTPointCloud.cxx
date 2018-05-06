#include "WireCellData/ToyCTPointCloud.h"

using namespace WireCell;

ToyCTPointCloud::ToyCTPointCloud(int u_min_ch, int u_max_ch, int v_min_ch, int v_max_ch, int w_min_ch, int w_max_ch, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u, double angle_v, double angle_w)
  : u_min_ch(u_min_ch)
  , u_max_ch(u_max_ch)
  , v_min_ch(v_min_ch)
  , v_max_ch(v_max_ch)
  , w_min_ch(w_min_ch)
  , w_max_ch(w_max_ch)
  , offset_t(offset_t)
  , offset_u(offset_u)
  , offset_v(offset_v)
  , offset_w(offset_w)
  , slope_t(slope_t)
  , slope_u(slope_u)
  , slope_v(slope_v)
  , slope_w(slope_w)
  , angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
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

void ToyCTPointCloud::AddPoint(int ch, int time_slice, int charge, int charge_err){
  int plane = 0;
  if (ch >= v_min_ch && ch <= v_max_ch){
    plane = 1;
  }else if (ch >= w_min_ch && ch <= w_max_ch){
    plane = 2;
  }
  
  CTPointCloud<double>::CTPoint point;

  point.channel = ch;
  point.time_slice = time_slice;
  point.charge = charge;
  point.charge_err = charge_err;

  point.x = (time_slice-offset_t)/slope_t;
  if (plane==0){
    point.y = (ch-offset_u)/slope_u;
    point.index = cloud_u.pts.size();
    cloud_u.pts.push_back(point);
  }else if (plane==1){
    point.y = (ch-offset_v)/slope_v;
    point.index = cloud_v.pts.size();
    cloud_v.pts.push_back(point);
  }else if (plane==2){
    point.y = (ch-offset_w)/slope_w;
    point.index = cloud_w.pts.size();
    cloud_w.pts.push_back(point);
  }
}


void ToyCTPointCloud::AddPoints(std::vector<int> *timesliceId, std::vector<std::vector<int>> *timesliceChannel, std::vector<std::vector<int>> *raw_charge , std::vector<std::vector<int>> *raw_charge_err){
  for (size_t i=0;i!=timesliceId->size(); i++){
    int time_slice = timesliceId->at(i);
    for (size_t j=0;j!=timesliceChannel->at(i).size();j++){
      int ch = timesliceChannel->at(i).at(j);
      int charge = raw_charge->at(i).at(j);
      int charge_err = raw_charge_err->at(i).at(j);
      AddPoint(time_slice,ch,charge,charge_err);
    }
  }
  
}

void ToyCTPointCloud::build_kdtree_index(){
  if (index_u!=0)    delete index_u;
  if (index_v!=0)    delete index_v;
  if (index_w!=0)    delete index_w;

  index_u = new my_kd_tree_ct_t(2, cloud_u, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_u->buildIndex();

  index_v = new my_kd_tree_ct_t(2, cloud_v, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_v->buildIndex();

  index_w = new my_kd_tree_ct_t(2, cloud_w, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_w->buildIndex();
}


int ToyCTPointCloud::get_num_points(int plane){
  if (plane==0){
    return cloud_u.pts.size();
  }else if (plane==1){
    return cloud_v.pts.size();
  }else if (plane==2){
    return cloud_w.pts.size();
  }
}

CTPointCloud<double>& ToyCTPointCloud::get_cloud(int plane){
  if (plane==0){
    return cloud_u;
  }else if (plane==1){
    return cloud_v;
  }else{
    return cloud_w;
  }
}

WireCell::CTPointCloud<double> ToyCTPointCloud::get_closest_points(WireCell::Point& p, double radius, int plane){
  WireCell::CTPointCloud<double> nearby_points;

  
  
  return nearby_points;
}
