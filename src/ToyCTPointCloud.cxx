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

void ToyCTPointCloud::Print(WireCell::Point &p){
  std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
  std::cout << p.x*slope_t+offset_t << std::endl;
  std::cout << (cos(angle_u) * p.z - sin(angle_u) *p.y)*slope_u+offset_u << std::endl;
  std::cout << (cos(angle_v) * p.z - sin(angle_v) *p.y)*slope_v+offset_v << std::endl;
  std::cout << (cos(angle_w) * p.z - sin(angle_w) *p.y)*slope_w+offset_w << std::endl;
  
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
    point.y = (ch-offset_u-u_min_ch)/slope_u;
    point.index = cloud_u.pts.size();
    cloud_u.pts.push_back(point);
  }else if (plane==1){
    point.y = (ch-offset_v-v_min_ch)/slope_v;
    point.index = cloud_v.pts.size();
    cloud_v.pts.push_back(point);
  }else if (plane==2){
    point.y = (ch-offset_w-w_min_ch)/slope_w;
    point.index = cloud_w.pts.size();
    cloud_w.pts.push_back(point);
  }
}

std::vector<int> ToyCTPointCloud::convert_3Dpoint_time_ch(WireCell::Point& p){
  std::vector<int> results(4,0);

  int time_slice = std::round(p.x * slope_t + offset_t);
  
  // U plane ... 
  double y = cos(angle_u) * p.z - sin(angle_u) *p.y;
  int ch_u = std::round(y * slope_u + offset_u + u_min_ch);
  
  y = cos(angle_v) * p.z - sin(angle_v) *p.y;
  int ch_v = std::round(y * slope_v + offset_v + v_min_ch);
  
  y = cos(angle_w) * p.z - sin(angle_w) *p.y;
  int ch_w = std::round(y*slope_w + offset_w + w_min_ch);

  if (ch_u < u_min_ch){
    ch_u = u_min_ch;
  }else if (ch_u > u_max_ch){
    ch_u = u_max_ch;
  }
  
  if (ch_v < v_min_ch){
    ch_v = v_min_ch;
  }else   if (ch_v > v_max_ch){
    ch_v = v_max_ch;
  }
  
  if (ch_w < w_min_ch) {
    ch_w = w_min_ch;
  }else if (ch_w > w_max_ch) {
    ch_w = w_max_ch;
  }
    
  results.at(0) = time_slice;
  results.at(1) = ch_u;
  results.at(2) = ch_v;
  results.at(3) = ch_w;
  
  return results;
}


void ToyCTPointCloud::AddPoints(std::vector<int> *timesliceId, std::vector<std::vector<int>> *timesliceChannel, std::vector<std::vector<int>> *raw_charge , std::vector<std::vector<int>> *raw_charge_err){
  for (size_t i=0;i!=timesliceId->size(); i++){
    int time_slice = timesliceId->at(i);
    for (size_t j=0;j!=timesliceChannel->at(i).size();j++){
      int ch = timesliceChannel->at(i).at(j);
      int charge = raw_charge->at(i).at(j);
      int charge_err = raw_charge_err->at(i).at(j);

      // std::cout << time_slice << " " << ch << " " << charge << " " << charge_err << std::endl;
      
      AddPoint(ch,time_slice,charge,charge_err);
    }
  }
  
}


void ToyCTPointCloud::AddDeadChs(std::map<int,std::pair<double,double> >& dead_u_index, std::map<int,std::pair<double,double> >& dead_v_index, std::map<int,std::pair<double,double> >& dead_w_index){
  dead_uchs = dead_u_index;
  dead_vchs = dead_v_index;
  dead_wchs = dead_w_index;
  
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


std::vector<std::pair<size_t,double>> ToyCTPointCloud::get_closest_index(WireCell::Point& p, double search_radius, int plane){
  double x, y;
  if (plane==0){
    x = p.x;
    y = cos(angle_u) * p.z - sin(angle_u) *p.y;
  }else if (plane==1){
    x = p.x;
    y = cos(angle_v) * p.z - sin(angle_v) *p.y;
  }else{
    x = p.x;
    y = cos(angle_w) * p.z - sin(angle_w) *p.y;
  }
  double query_pt[2];
  query_pt[0] = x;
  query_pt[1] = y;

  std::vector<std::pair<size_t,double> >   ret_matches;
  nanoflann::SearchParams params;
  if (plane ==0){
    const size_t nMatches = index_u->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }else if (plane==1){
    const size_t nMatches = index_v->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }else{
    const size_t nMatches = index_w->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }
  return ret_matches;
}

bool ToyCTPointCloud::is_good_point(WireCell::Point& p, double radius, int ch_range, int allowed_bad){

 
  
  int num_planes = 0;
  {
    WireCell::CTPointCloud<double> pts = get_closest_points(p, radius, 0);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 0, ch_range))
	num_planes ++;
    }
  }
  {
    WireCell::CTPointCloud<double> pts = get_closest_points(p, radius, 1);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 1, ch_range))
	num_planes ++;
    }
  }
  {
    WireCell::CTPointCloud<double> pts = get_closest_points(p, radius, 2);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 2, ch_range))
	num_planes ++;
    }
  }
  if (num_planes >= 3 - allowed_bad)
    return true;

  return false;
  
}


bool ToyCTPointCloud::get_closest_dead_chs(WireCell::Point& p, int plane, int ch_range){

  std::vector<int> results = convert_3Dpoint_time_ch(p);

  results.at(2) -= u_max_ch - u_min_ch + 1;
  results.at(3) -= v_max_ch - v_min_ch + 1 + u_max_ch - u_min_ch + 1;
  
  // if (p.x > 275*units::cm && p.x < 290*units::cm && p.y > 25*units::cm &&
  //     p.y < 35*units::cm && p.z > 210*units::cm && p.z < 230*units::cm){
  //   bool flag_u = dead_uchs.find(results.at(1))==dead_uchs.end();
  //   bool flag_v = dead_vchs.find(results.at(2))==dead_vchs.end();
  //   bool flag_w = dead_wchs.find(results.at(3))==dead_wchs.end();
  //   std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " << results.at(0) << " " << results.at(1) << " " << results.at(2) << " " << " " << results.at(3) << " " << flag_u << " " << flag_v  << " " << flag_w << " " << std::endl;
  //   if (!flag_u)
  //     std::cout << "U " << dead_uchs[results.at(1)].first/units::cm << " " << dead_uchs[results.at(1)].second/units::cm << std::endl;
  //   if (!flag_v)
  //     std::cout << "V " << dead_vchs[results.at(2)].first/units::cm << " " << dead_vchs[results.at(2)].second/units::cm << std::endl;
    
  //   // for (auto it = dead_vchs.begin(); it!=dead_vchs.end(); it++){
  //   //   std::cout << it->first << " " << it->second.first/units::cm << " " << it->second.second/units::cm << std::endl;
  //   // }
  // }

  
  
  if (plane == 0){
    for (int ch = results.at(1) - ch_range; ch <= results.at(1) + ch_range; ch ++){
      if (dead_uchs.find(ch)==dead_uchs.end()) continue;
      if (p.x >= dead_uchs[ch].first && p.x <= dead_uchs[ch].second)
	return true;
    }
  }else if (plane == 1){
    for (int ch = results.at(2) - ch_range; ch <= results.at(2) + ch_range; ch ++){
      if (dead_vchs.find(ch)==dead_vchs.end()) continue;
      if (p.x >= dead_vchs[ch].first && p.x <= dead_vchs[ch].second)
	return true;
    }
  }else if (plane == 2){
    for (int ch = results.at(3) - ch_range; ch <= results.at(3) + ch_range; ch ++){
      if (dead_wchs.find(ch)==dead_wchs.end()) continue;
      if (p.x >= dead_wchs[ch].first && p.x <= dead_wchs[ch].second)
	return true;
    }
  }

  
  return false;
}

WireCell::CTPointCloud<double> ToyCTPointCloud::get_closest_points(WireCell::Point& p, double radius, int plane){
  WireCell::CTPointCloud<double> nearby_points;

  std::vector<std::pair<size_t,double>> indices = get_closest_index(p, radius, plane);
  for (size_t i=0;i!=indices.size();i++){
    if (plane==0){
      nearby_points.pts.push_back(cloud_u.pts.at(indices.at(i).first));
    }else if (plane==1){
      nearby_points.pts.push_back(cloud_v.pts.at(indices.at(i).first));
    }else if (plane==2){
      nearby_points.pts.push_back(cloud_w.pts.at(indices.at(i).first));
    }
  }
  
  
  return nearby_points;
}


