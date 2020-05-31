#include "WCPData/ToyCTPointCloud.h"

using namespace WCP;

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

void ToyCTPointCloud::Print(WCP::Point &p){
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

std::pair<double,double> ToyCTPointCloud::convert_time_ch_2Dpoint(int timeslice, int channel, int plane){
  double x = (timeslice - offset_t)/slope_t;
  double y;
  if (plane==0){
    y = (channel - u_min_ch - offset_u)/slope_u;
  }else if(plane==1){
    y = (channel - v_min_ch - offset_v)/slope_v;
  }else if(plane==2){
    y = (channel - w_min_ch - offset_w)/slope_w;
  }
  return std::make_pair(x,y);
}

std::vector<int> ToyCTPointCloud::convert_3Dpoint_time_ch(WCP::Point& p){
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

void ToyCTPointCloud::UpdateDeadChs(){
  std::map<int, std::pair<double, double> > live_uchs;
  std::map<int, std::pair<double, double> > live_vchs;
  std::map<int, std::pair<double, double> > live_wchs;
  for (size_t i=0;i!=cloud_u.pts.size();i++){
    int ch = std::round(cloud_u.pts.at(i).channel);
    double xpos = (cloud_u.pts.at(i).time_slice - offset_t)/slope_t;
    if (live_uchs.find(ch) == live_uchs.find(ch)){
      live_uchs[ch] = std::make_pair(xpos - 0.1*units::cm, xpos + 0.1*units::cm);
    }else{
      if (xpos + 0.1*units::cm > live_uchs[ch].second) live_uchs[ch].second = xpos+0.1*units::cm;
      if (xpos - 0.1*units::cm < live_uchs[ch].first) live_uchs[ch].first = xpos - 0.1*units::cm;
    }
  }
  for (size_t i=0;i!=cloud_v.pts.size();i++){
    int ch = std::round(cloud_v.pts.at(i).channel)-2400;
    double xpos = (cloud_v.pts.at(i).time_slice - offset_t)/slope_t;
    if (live_vchs.find(ch) == live_vchs.find(ch)){
      live_vchs[ch] = std::make_pair(xpos - 0.1*units::cm, xpos + 0.1*units::cm);
    }else{
      if (xpos + 0.1*units::cm > live_vchs[ch].second) live_vchs[ch].second = xpos+0.1*units::cm;
      if (xpos - 0.1*units::cm < live_vchs[ch].first) live_vchs[ch].first = xpos - 0.1*units::cm;
    }
  }
 for (size_t i=0;i!=cloud_w.pts.size();i++){
    int ch = std::round(cloud_w.pts.at(i).channel)-4800;
    double xpos = (cloud_w.pts.at(i).time_slice - offset_t)/slope_t;
    if (live_wchs.find(ch) == live_wchs.find(ch)){
      live_wchs[ch] = std::make_pair(xpos - 0.1*units::cm, xpos + 0.1*units::cm);
    }else{
      if (xpos + 0.1*units::cm > live_wchs[ch].second) live_wchs[ch].second = xpos+0.1*units::cm;
      if (xpos - 0.1*units::cm < live_wchs[ch].first) live_wchs[ch].first = xpos - 0.1*units::cm;
    }
  }

 for (auto it = live_vchs.begin(); it != live_vchs.end(); it++){
   int ch = it->first;
   //   std::cout << "A: " << ch << std::endl;

   if (dead_vchs.find(ch) != dead_vchs.end())
     std::cout << ch << " " << it->second.first << " " << it->second.second << " " << dead_vchs[ch].first << " " << dead_vchs[ch].second << std::endl;
 }
 

 for (auto it = dead_vchs.begin(); it != dead_vchs.end(); it++){
   int ch = it->first;
   std::cout << "B: " << ch << std::endl;
 }
  
}

void ToyCTPointCloud::AddDeadChs(std::map<int,std::pair<double,double> >& dead_u_index, std::map<int,std::pair<double,double> >& dead_v_index, std::map<int,std::pair<double,double> >& dead_w_index){
  dead_uchs = dead_u_index;
  dead_vchs = dead_v_index;
  dead_wchs = dead_w_index;


  // examine a bit dead channels ...
  //double min_xpos = (min_time-offset_t)/slope_t;
  //double max_xpos = (max_time-offset_t)/slope_t;
  // point.channel = ch;
  //point.time_slice = time_slice;
  
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


std::vector<std::pair<size_t,double>> ToyCTPointCloud::get_closest_index(WCP::Point& p, double search_radius, int plane){
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
  //  std::cout << x << " " << y << " " << plane << std::endl;
  
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

double ToyCTPointCloud::get_ave_3d_charge(WCP::Point& p, double radius){
  double charge = 0;
  int ncount = 0;
  if (!get_closest_dead_chs(p,0)){
    charge += get_ave_charge(p, radius,0);
    ncount ++;
  }
  if (!get_closest_dead_chs(p,1)){
    charge += get_ave_charge(p, radius,1);
    ncount ++;
  }
  if (!get_closest_dead_chs(p,2)){
    charge += get_ave_charge(p, radius,2);
    ncount ++;
  }
  if (ncount!=0) charge /= ncount;
  return charge;
}


double ToyCTPointCloud::get_ave_charge(WCP::Point& p, double radius, int plane){
  double sum_charge = 0;
  double ncount = 0;
  WCP::CTPointCloud<double> pts = get_closest_points(p, radius, plane);
  // std::cout << pts.pts.size() << std::endl;
  for (size_t i=0;i!=pts.pts.size(); i++){
    sum_charge += pts.pts.at(i).charge;
    ncount ++;
  }
  if (ncount!=0) sum_charge/=ncount;
  return sum_charge;
}

std::vector<int>  ToyCTPointCloud::test_good_point(WCP::Point& p, double radius, int ch_range){
  std::vector<int> num_planes;//[6]={0,0,0,0,0,0};
  num_planes.resize(6,0);
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 0);
    if (pts.pts.size()>0){
      num_planes[0] ++;
    }else{
      if (get_closest_dead_chs(p, 0, ch_range))
	num_planes[3] ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 1);
    if (pts.pts.size()>0){
      num_planes[1] ++;
    }else{
      if (get_closest_dead_chs(p, 1, ch_range))
	num_planes[4] ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 2);
    if (pts.pts.size()>0){
      num_planes[2] ++;
    }else{
      if (get_closest_dead_chs(p, 2, ch_range))
	num_planes[5] ++;
    }
  }

  return num_planes;//std::make_tuple(num_planes[0],num_planes[1],num_planes[2]);
}

bool ToyCTPointCloud::is_good_point_wc(WCP::Point& p, double radius, int ch_range, int allowed_bad){
  int num_planes = 0;
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 0);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 0, ch_range))
	num_planes ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 1);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 1, ch_range))
	num_planes ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 2);
    if (pts.pts.size()>0){
      num_planes +=2;
    }else{
      if (get_closest_dead_chs(p, 2, ch_range))
	num_planes +=2;
    }
  }
  if (num_planes >= 4 - allowed_bad)
    return true;

  return false;
  
}

bool ToyCTPointCloud::is_good_point(WCP::Point& p, double radius, int ch_range, int allowed_bad){
  int num_planes = 0;
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 0);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 0, ch_range))
	num_planes ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 1);
    if (pts.pts.size()>0){
      num_planes ++;
    }else{
      if (get_closest_dead_chs(p, 1, ch_range))
	num_planes ++;
    }
  }
  {
    WCP::CTPointCloud<double> pts = get_closest_points(p, radius, 2);
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


std::map<std::pair<int,int>, std::pair<double,double> > ToyCTPointCloud::get_overlap_good_ch_charge(int min_time, int max_time, int min_ch, int max_ch, int plane_no){
  int ch_offset = 0;
  if (plane_no == 1){
    ch_offset = u_max_ch - u_min_ch + 1;
  }else if (plane_no==2){
    ch_offset = v_max_ch - v_min_ch + 1 + u_max_ch - u_min_ch + 1;
  }
  // time channel, charge, charge_err ...
  std::map<std::pair<int,int>, std::pair<double, double> > map_time_ch_charge;
  if (plane_no==0){
    for (size_t i=0;i!=cloud_u.pts.size();i++){
       
      if (cloud_u.pts.at(i).time_slice>=min_time && cloud_u.pts.at(i).time_slice <= max_time &&
	  cloud_u.pts.at(i).channel >= min_ch && cloud_u.pts.at(i).channel <= max_ch)
	map_time_ch_charge[std::make_pair(cloud_u.pts.at(i).time_slice, cloud_u.pts.at(i).channel)] = std::make_pair(cloud_u.pts.at(i).charge, cloud_u.pts.at(i).charge_err);
    }
	   
  }else if (plane_no == 1){
    for (size_t i=0;i!=cloud_v.pts.size();i++){
      if (cloud_v.pts.at(i).time_slice>=min_time && cloud_v.pts.at(i).time_slice <= max_time &&
	  cloud_v.pts.at(i).channel >= min_ch && cloud_v.pts.at(i).channel <= max_ch)
	map_time_ch_charge[std::make_pair(cloud_v.pts.at(i).time_slice, cloud_v.pts.at(i).channel)] = std::make_pair(cloud_v.pts.at(i).charge, cloud_v.pts.at(i).charge_err);
    }
  }else if (plane_no == 2){
    for (size_t i=0;i!=cloud_w.pts.size();i++){
      //std::cout << cloud_w.pts.at(i).channel << " haha " << min_ch << " " << max_ch << std::endl;
      if (cloud_w.pts.at(i).time_slice>=min_time && cloud_w.pts.at(i).time_slice <= max_time &&
	  cloud_w.pts.at(i).channel >= min_ch && cloud_w.pts.at(i).channel <= max_ch)
	map_time_ch_charge[std::make_pair(cloud_w.pts.at(i).time_slice, cloud_w.pts.at(i).channel)] = std::make_pair(cloud_w.pts.at(i).charge, cloud_w.pts.at(i).charge_err);
    }
  }
  
  return map_time_ch_charge;
  
}

std::map<int, std::pair<int, int> > ToyCTPointCloud::get_all_dead_chs(){
  std::map<int, std::pair<int, int> > results;
  
  int ch_offset = 0;
  for (auto it = dead_uchs.begin(); it!=dead_uchs.end(); it++){
    int temp_ch = it->first+ch_offset;
    int min_time = it->second.first * slope_t + offset_t - 3;
    int max_time = it->second.second * slope_t + offset_t + 3;
    results[temp_ch] = std::make_pair(min_time, max_time);
  }
  ch_offset = u_max_ch - u_min_ch + 1;
  for (auto it = dead_vchs.begin(); it!=dead_vchs.end(); it++){
    int temp_ch = it->first+ch_offset;
    int min_time = it->second.first * slope_t + offset_t - 3;
    int max_time = it->second.second * slope_t + offset_t + 3;
    results[temp_ch] = std::make_pair(min_time, max_time);
  }
  ch_offset = v_max_ch - v_min_ch + 1 + u_max_ch - u_min_ch + 1;
  for (auto it = dead_wchs.begin(); it!=dead_wchs.end(); it++){
    int temp_ch = it->first+ch_offset;
    int min_time = it->second.first * slope_t + offset_t - 3;
    int max_time = it->second.second * slope_t + offset_t + 3;
    results[temp_ch] = std::make_pair(min_time, max_time);
  }

  return results;
}

std::vector<std::pair<int, int> > ToyCTPointCloud::get_overlap_dead_chs(int min_time, int max_time, int min_ch, int max_ch, int plane_no, bool flag_ignore_time){
  int ch_offset = 0;
  if (plane_no == 1){
    ch_offset = u_max_ch - u_min_ch + 1;
  }else if (plane_no==2){
    ch_offset = v_max_ch - v_min_ch + 1 + u_max_ch - u_min_ch + 1;
  }

  double min_xpos = (min_time-offset_t)/slope_t;
  double max_xpos = (max_time-offset_t)/slope_t;
  
  std::set<int> dead_chs;
  // U plane
  if (plane_no==0){
    for (auto it = dead_uchs.begin(); it!=dead_uchs.end(); it++){
      int temp_ch = it->first+ch_offset;
      double temp_min_xpos = it->second.first;
      double temp_max_xpos = it->second.second;
      if (flag_ignore_time){
	if (temp_ch>=min_ch && temp_ch<=max_ch)
	  dead_chs.insert(temp_ch);
      }else{
	if (temp_ch>=min_ch && temp_ch<=max_ch &&
	    max_xpos >= temp_min_xpos && min_xpos <= temp_max_xpos)
	  dead_chs.insert(temp_ch);
      }
    }
    
  }else if (plane_no==1){ // V plane
    for (auto it = dead_vchs.begin(); it!=dead_vchs.end(); it++){
      int temp_ch = it->first+ch_offset;
      double temp_min_xpos = it->second.first;
      double temp_max_xpos = it->second.second;
      if (flag_ignore_time){
	if (temp_ch>=min_ch && temp_ch<=max_ch)
	  dead_chs.insert(temp_ch);
      }else{
	if (temp_ch>=min_ch && temp_ch<=max_ch &&
	    max_xpos >= temp_min_xpos && min_xpos <= temp_max_xpos)
	  dead_chs.insert(temp_ch);
      }
    }
    
  }else if (plane_no==2){ // W plane 
    for (auto it = dead_wchs.begin(); it!=dead_wchs.end(); it++){
      int temp_ch = it->first+ch_offset;
      double temp_min_xpos = it->second.first;
      double temp_max_xpos = it->second.second;
      //  std::cout << temp_ch << " " << min_ch << " " << max_ch << std::endl;
      if (flag_ignore_time){
	if (temp_ch>=min_ch && temp_ch<=max_ch)
	  dead_chs.insert(temp_ch);
      }else{
	if (temp_ch>=min_ch && temp_ch<=max_ch &&
	    max_xpos >= temp_min_xpos && min_xpos <= temp_max_xpos)
	  dead_chs.insert(temp_ch);
      }
    }
    
  }

  std::vector<std::pair<int, int> > dead_ch_range;
  for (auto it = dead_chs.begin(); it!=dead_chs.end(); it++){
    if (dead_ch_range.size()==0){
      dead_ch_range.push_back(std::make_pair(*it, *it));
    }else{
      if (*it-dead_ch_range.back().second==1){
	dead_ch_range.back().second = *it;
      }else{
	dead_ch_range.push_back(std::make_pair(*it, *it));
      }
    }
  }
  
  return dead_ch_range;
}




bool ToyCTPointCloud::get_closest_dead_chs(WCP::Point& p, int plane, int ch_range){

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

WCP::CTPointCloud<double> ToyCTPointCloud::get_closest_points(WCP::Point& p, double radius, int plane){
  WCP::CTPointCloud<double> nearby_points;

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


