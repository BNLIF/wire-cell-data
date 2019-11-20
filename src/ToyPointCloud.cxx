#include "WCPData/ToyPointCloud.h"

using namespace WCP;

WCP::ToyPointCloud::ToyPointCloud(double angle_u, double angle_v, double angle_w)
  : angle_u(angle_u)
  , angle_v(angle_v)
  , angle_w(angle_w)
{
  index = 0;
}

WCP::ToyPointCloud::~ToyPointCloud(){
  if (index!=0){
    delete index;
    delete index_u;
    delete index_v;
    delete index_w;
  }
  map_mcell_indices.clear();
  cloud.pts.clear();
  cloud_u.pts.clear();
  cloud_v.pts.clear();
  cloud_w.pts.clear();
}

std::pair<int,double> WCP::ToyPointCloud::get_closest_point_along_vec(Point& p_test1, TVector3& dir, double test_dis, double dis_step, double angle_cut, double dis_cut){
  if (dir.Mag()!=1)
    dir.SetMag(1);
  
  bool flag = false;
  Point p_test;
  double min_dis = 1e9;
  double min_dis1 = 1e9;
  int min_index = -1;
  for (int i=0; i!= int(test_dis/dis_step)+1;i++){
    p_test.x = p_test1.x + dir.X() * i * dis_step;
    p_test.y = p_test1.y + dir.Y() * i * dis_step;
    p_test.z = p_test1.z + dir.Z() * i * dis_step;
    
    WCP::WCPointCloud<double>::WCPoint& closest_pt = get_closest_wcpoint(p_test);
    double dis = sqrt(pow(p_test.x - closest_pt.x,2)+pow(p_test.y - closest_pt.y,2)+pow(p_test.z - closest_pt.z,2));
    double dis1 = sqrt(pow(p_test1.x - closest_pt.x,2)+pow(p_test1.y - closest_pt.y,2)+pow(p_test1.z - closest_pt.z,2));

    // if (fabs(p_test1.x-199.942*units::cm) < 1*units::cm &&
    // 	fabs(p_test1.y+19.4682*units::cm) < 1*units::cm &&
    // 	fabs(p_test1.z-127*units::cm) < 1*units::cm)
    //   std::cout << i << " " << dis/units::cm << " " << std::min(dis1 * tan(angle_cut/180.*3.1415926),dis_cut)/units::cm << " " << p_test.x/units::cm << " " << p_test.y/units::cm << " " << p_test.z/units::cm << std::endl;
	      
    
    if (dis < std::min(dis1 * tan(angle_cut/180.*3.1415926),dis_cut)){
      if (dis < min_dis){
	min_dis = dis;
	min_index = closest_pt.index;
	min_dis1 = dis1;
      }
      if (dis < 3*units::cm)
	return std::make_pair(closest_pt.index,dis1);
    }
  }

  return std::make_pair(min_index,min_dis1);
}

std::tuple<int,int,double> WCP::ToyPointCloud::get_closest_points(ToyPointCloud *point_cloud){

  //std::cout << cloud.pts.size() << " " << point_cloud->get_cloud().pts.size() << std::endl;
   
  
  WCP::WCPointCloud<double>::WCPoint p1 = cloud.pts.front();
  WCP::WCPointCloud<double>::WCPoint p2 = point_cloud->get_cloud().pts.front();
  WCP::WCPointCloud<double>::WCPoint p1_save;
  WCP::WCPointCloud<double>::WCPoint p2_save;
  double min_dis = 1e9;
  
  int prev_index1 = -1;
  int prev_index2 = -1;
  while(p1.index!=prev_index1 || p2.index!=prev_index2){
    prev_index1 = p1.index;
    prev_index2 = p2.index;
    p2 = point_cloud->get_closest_wcpoint(p1);
    p1 = get_closest_wcpoint(p2);
  }
  double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  if (dis < min_dis){
    min_dis = dis;
    p1_save = p1;
    p2_save = p2;
  }

  prev_index1 = -1;
  prev_index2 = -1;
  p1 = cloud.pts.back();
  p2 = point_cloud->get_cloud().pts.front();
  while(p1.index!=prev_index1 || p2.index!=prev_index2){
    prev_index1 = p1.index;
    prev_index2 = p2.index;
    p2 = point_cloud->get_closest_wcpoint(p1);
    p1 = get_closest_wcpoint(p2);
  }
  dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  if (dis < min_dis){
    min_dis = dis;
    p1_save = p1;
    p2_save = p2;
  }

  prev_index1 = -1;
  prev_index2 = -1;
  p1 = cloud.pts.back();
  p2 = point_cloud->get_cloud().pts.front();
  while(p1.index!=prev_index1 || p2.index!=prev_index2){
    prev_index1 = p1.index;
    prev_index2 = p2.index;
    p1 = get_closest_wcpoint(p2);
    p2 = point_cloud->get_closest_wcpoint(p1);
  }
  dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  if (dis < min_dis){
    min_dis = dis;
    p1_save = p1;
    p2_save = p2;
  }

  prev_index1 = -1;
  prev_index2 = -1;
  p1 = cloud.pts.back();
  p2 = point_cloud->get_cloud().pts.back();
  while(p1.index!=prev_index1 || p2.index!=prev_index2){
    prev_index1 = p1.index;
    prev_index2 = p2.index;
    p1 = get_closest_wcpoint(p2);
    p2 = point_cloud->get_closest_wcpoint(p1);
  }
  dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
  if (dis < min_dis){
    min_dis = dis;
    p1_save = p1;
    p2_save = p2;
  }

  
  
  
  return std::make_tuple(p1_save.index,p2_save.index,min_dis);
}

WCP::WCPointCloud<double>::WCPoint& WCP::ToyPointCloud::get_closest_wcpoint(WCP::WCPointCloud<double>::WCPoint& wcp){
  Point p(wcp.x,wcp.y,wcp.z);
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,1);
  return cloud.pts[results.front().first];
}

WCP::WCPointCloud<double>::WCPoint& WCP::ToyPointCloud::get_closest_wcpoint(WCP::Point& p){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,1);
  return cloud.pts[results.front().first];
}


void WCP::ToyPointCloud::AddPoint(WCPointCloud<double>::WCPoint& wcp, WC2DPointCloud<double>::WC2DPoint& wcp_u, WC2DPointCloud<double>::WC2DPoint& wcp_v, WC2DPointCloud<double>::WC2DPoint& wcp_w){
  SlimMergeGeomCell *mcell = wcp.mcell;
  int index = wcp.index;
  
  if (map_mcell_indices.find(mcell)==map_mcell_indices.end()){
    std::vector<int> temp_indices;
    temp_indices.push_back(index);
    map_mcell_indices[mcell] = temp_indices;
  }else{
    map_mcell_indices[mcell].push_back(index);
  }
  

  cloud.pts.push_back(wcp);
  cloud_u.pts.push_back(wcp_u);
  cloud_v.pts.push_back(wcp_v);
  cloud_w.pts.push_back(wcp_w);

  
  
}


void WCP::ToyPointCloud::AddPoint(Point& p, std::tuple<int,int,int>& wires_index, SlimMergeGeomCell *mcell){
  WCPointCloud<double>::WCPoint point;
  point.x = p.x;
  point.y = p.y;
  point.z = p.z;

  point.index_u = std::get<0>(wires_index);
  point.index_v = std::get<1>(wires_index);
  point.index_w = std::get<2>(wires_index);
    
  point.mcell = mcell;
  point.index =cloud.pts.size();

  WC2DPointCloud<double>::WC2DPoint point_u;
  point_u.x = p.x;
  point_u.y = cos(angle_u) * p.z - sin(angle_u) *p.y;
  point_u.index = cloud_u.pts.size();
  
  WC2DPointCloud<double>::WC2DPoint point_v;
  point_v.x = p.x;
  point_v.y = cos(angle_v) * p.z - sin(angle_v) *p.y;
  point_v.index = cloud_v.pts.size();
  
  WC2DPointCloud<double>::WC2DPoint point_w;
  point_w.x = p.x;
  point_w.y = cos(angle_w) * p.z - sin(angle_w) *p.y;
  point_w.index = cloud_w.pts.size();
  
  if (map_mcell_indices.find(mcell)==map_mcell_indices.end()){
    std::vector<int> temp_indices;
    temp_indices.push_back(point.index);
    map_mcell_indices[mcell] = temp_indices;
  }else{
    map_mcell_indices[mcell].push_back(point.index);
  }
  

  cloud.pts.push_back(point);
  cloud_u.pts.push_back(point_u);
  cloud_v.pts.push_back(point_v);
  cloud_w.pts.push_back(point_w);
  
}



void WCP::ToyPointCloud::AddPoints(PointVector& ps){
  size_t current_size = cloud.pts.size();
  cloud.pts.resize(current_size + ps.size());
  cloud_u.pts.resize(current_size + ps.size());
  cloud_v.pts.resize(current_size + ps.size());
  cloud_w.pts.resize(current_size + ps.size());
  
  for (size_t i=0;i!=ps.size();i++){
    cloud.pts[current_size+i].x = ps.at(i).x;
    cloud.pts[current_size+i].y = ps.at(i).y;
    cloud.pts[current_size+i].z = ps.at(i).z;
    cloud.pts[current_size+i].index_u = 0;
    cloud.pts[current_size+i].index_v = 0;
    cloud.pts[current_size+i].index_w = 0;
    cloud.pts[current_size+i].mcell = 0;
    cloud.pts[current_size+i].index = current_size+i;

    cloud_u.pts[current_size+i].x = ps.at(i).x;
    cloud_u.pts[current_size+i].y = cos(angle_u) * ps.at(i).z - sin(angle_u) *ps.at(i).y;
    cloud_u.pts[current_size+i].index = current_size+i;
    
    cloud_v.pts[current_size+i].x = ps.at(i).x;
    cloud_v.pts[current_size+i].y = cos(angle_v) * ps.at(i).z - sin(angle_v) *ps.at(i).y;
    cloud_v.pts[current_size+i].index = current_size+i;
    
    cloud_w.pts[current_size+i].x = ps.at(i).x;
    cloud_w.pts[current_size+i].y = cos(angle_w) * ps.at(i).z - sin(angle_w) *ps.at(i).y;
    cloud_w.pts[current_size+i].index = current_size+i;
    

  }
}


void WCP::ToyPointCloud::AddPoints(PointVector& ps, std::vector<std::tuple<int,int,int>>& wires_indices, SlimMergeGeomCell *mcell){
  size_t current_size = cloud.pts.size();
  cloud.pts.resize(current_size + ps.size());
  cloud_u.pts.resize(current_size + ps.size());
  cloud_v.pts.resize(current_size + ps.size());
  cloud_w.pts.resize(current_size + ps.size());
  
  
  if (map_mcell_indices.find(mcell)==map_mcell_indices.end()){
    std::vector<int> temp_pts;
    map_mcell_indices[mcell] = temp_pts;
  }
  
  for (size_t i=0;i!=ps.size();i++){
    cloud.pts[current_size+i].x = ps.at(i).x;
    cloud.pts[current_size+i].y = ps.at(i).y;
    cloud.pts[current_size+i].z = ps.at(i).z;
    cloud.pts[current_size+i].index_u = std::get<0>(wires_indices.at(i)) ;
    cloud.pts[current_size+i].index_v = std::get<1>(wires_indices.at(i)) ;
    cloud.pts[current_size+i].index_w = std::get<2>(wires_indices.at(i)) ;
    cloud.pts[current_size+i].mcell = mcell;
    cloud.pts[current_size+i].index = current_size+i;

    cloud_u.pts[current_size+i].x = ps.at(i).x;
    cloud_u.pts[current_size+i].y = cos(angle_u) * ps.at(i).z - sin(angle_u) *ps.at(i).y;
    cloud_u.pts[current_size+i].index = current_size+i;
    
    cloud_v.pts[current_size+i].x = ps.at(i).x;
    cloud_v.pts[current_size+i].y = cos(angle_v) * ps.at(i).z - sin(angle_v) *ps.at(i).y;
    cloud_v.pts[current_size+i].index = current_size+i;
    
    cloud_w.pts[current_size+i].x = ps.at(i).x;
    cloud_w.pts[current_size+i].y = cos(angle_w) * ps.at(i).z - sin(angle_w) *ps.at(i).y;
    cloud_w.pts[current_size+i].index = current_size+i;
    

    map_mcell_indices[mcell].push_back(current_size+i);
  }
}

void WCP::ToyPointCloud::Print(){
  for (size_t i=0; i!=cloud.pts.size(); i++){
    std::cout << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " << cloud.pts[i].mcell << std::endl;
  }
}

void WCP::ToyPointCloud::build_kdtree_index(){
  if (index!=(my_kd_tree_t*)0){
    delete index;
    delete index_u;
    delete index_v;
    delete index_w;
  }
  index = new my_kd_tree_t(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();

  index_u = new my_kd_tree_2d_t(2, cloud_u, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_u->buildIndex();

  index_v = new my_kd_tree_2d_t(2, cloud_v, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_v->buildIndex();

  index_w = new my_kd_tree_2d_t(2, cloud_w, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index_w->buildIndex();
  
}



std::vector<std::pair<size_t,double>> WCP::ToyPointCloud::get_closest_index(WCP::Point& p, int N){
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;  
  
  N = index->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  ret_index.resize(N);
  out_dist_sqr.resize(N);
  
  std::vector<std::pair<size_t,double>> results(N);
  for(size_t i=0; i!=N; i++){
    results.at(i) = std::make_pair(ret_index.at(i),out_dist_sqr.at(i));
  }
  
  return results;
}



std::vector<std::pair<size_t,double>> WCP::ToyPointCloud::get_closest_2d_index(double x, double y, int N, int plane){
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  
  double query_pt[2];
  query_pt[0] = x;
  query_pt[1] = y;
  if (plane==0){
    N = index_u->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else if (plane==1){
    N = index_v->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else{
    N = index_w->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }
  ret_index.resize(N);
  out_dist_sqr.resize(N);
  
  std::vector<std::pair<size_t,double>> results(N);
  for(size_t i=0; i!=N; i++){
    results.at(i) = std::make_pair(ret_index.at(i),out_dist_sqr.at(i));
  }
  
  return results;
}

double WCP::ToyPointCloud::get_closest_2d_dis(double x, double y, int plane){
  double query_pt[2];
  query_pt[0] = x;
  query_pt[1] = y;
  int N = 1;
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  if (plane==0){
    N = index_u->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else if (plane==1){
    N = index_v->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else{
    N = index_w->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }
  ret_index.resize(N);
  out_dist_sqr.resize(N);

  if (N>0)
    return sqrt(out_dist_sqr[0]);
  else
    return 1e9;
}

std::pair<int,double> WCP::ToyPointCloud::get_closest_2d_dis(WCP::Point& p, int plane){
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
  int N = 1;
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  if (plane==0){
    N = index_u->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else if (plane==1){
    N = index_v->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else{
    N = index_w->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }
  ret_index.resize(N);
  out_dist_sqr.resize(N);

  if (N>0)
    return std::make_pair(ret_index[0],sqrt(out_dist_sqr[0]));
  else
    return std::make_pair(-1,1e9);
}


std::vector<size_t> WCP::ToyPointCloud::get_closest_2d_index(Point& p, int N, int plane){
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

  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  if (plane==0){
    N = index_u->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else if (plane==1){
    N = index_v->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }else{
    N = index_w->knnSearch(&query_pt[0], N, &ret_index[0], &out_dist_sqr[0]);
  }
  ret_index.resize(N);
  out_dist_sqr.resize(N);
  
  return ret_index;
}




std::vector<std::pair<size_t,double>> WCP::ToyPointCloud::get_closest_index(WCP::Point& p, double search_radius){
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;
  std::vector<std::pair<size_t,double> >   ret_matches;
  nanoflann::SearchParams params;
  const size_t nMatches = index->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  
  return ret_matches;
}


std::vector<std::pair<size_t,double>> WCP::ToyPointCloud::get_closest_2d_index(double x, double y, double search_radius, int plane){
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

std::vector<size_t> WCP::ToyPointCloud::get_closest_2d_index(Point& p, double search_radius, int plane){
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
  std::vector<std::pair<size_t,double>>   ret_matches;
  nanoflann::SearchParams params;
  if (plane ==0){
    const size_t nMatches = index_u->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }else if (plane==1){
    const size_t nMatches = index_v->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }else{
    const size_t nMatches = index_w->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  }
  std::vector<size_t> ret_index(ret_matches.size());
  for (size_t i=0;i!=ret_matches.size();i++){
    ret_index.at(i) = ret_matches.at(i).first;
  }
  
  return ret_index;
}


std::map<WCP::SlimMergeGeomCell*, Point> WCP::ToyPointCloud::get_closest_mcell(WCP::Point& p, int N){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,N);
  
  std::map<WCP::SlimMergeGeomCell*, Point> mcell_point_map;
  std::map<WCP::SlimMergeGeomCell*, double> mcell_dis_map;
  
  for (auto it = results.begin(); it!= results.end(); it++){
    size_t index = (*it).first;
    double dis = (*it).second;
    SlimMergeGeomCell *mcell = cloud.pts[index].mcell;
    Point p;
    p.x = cloud.pts[index].x;
    p.y = cloud.pts[index].y;
    p.z = cloud.pts[index].z;
    
    if (mcell_dis_map.find(mcell)==mcell_dis_map.end()){
      mcell_dis_map[mcell]=dis;
      mcell_point_map[mcell] = p;
    }else{
      if (dis < mcell_dis_map[mcell]){
	mcell_dis_map[mcell]=dis;
	mcell_point_map[mcell]=p;
      }
    }
  }
  
  return mcell_point_map;
}


std::map<WCP::SlimMergeGeomCell*, Point> WCP::ToyPointCloud::get_closest_mcell(WCP::Point& p, double search_radius){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,search_radius);
  
  std::map<WCP::SlimMergeGeomCell*, Point> mcell_point_map;
  std::map<WCP::SlimMergeGeomCell*, double> mcell_dis_map;
  
  for (auto it = results.begin(); it!= results.end(); it++){
    size_t index = (*it).first;
    double dis = (*it).second;

    //    std::cout << index << " " << dis/units::cm << std::endl;
    
    SlimMergeGeomCell *mcell = cloud.pts[index].mcell;
    Point p1;
    p1.x = cloud.pts[index].x;
    p1.y = cloud.pts[index].y;
    p1.z = cloud.pts[index].z;

    // Point pc = mcell->center();
    // std::cout << index << " " << dis/units::cm << " " << sqrt(pow(p1.x-p.x,2)+pow(p1.y-p.y,2)+pow(p1.z-p.z,2))/units::cm << " " << sqrt(pow(pc.x-p.x,2)+pow(pc.y-p.y,2)+pow(pc.z-p.z,2))/units::cm << std::endl;
    
    if (mcell_dis_map.find(mcell)==mcell_dis_map.end()){
      mcell_dis_map[mcell]=dis;
      mcell_point_map[mcell] = p1;
    }else{
      if (dis < mcell_dis_map[mcell]){
	mcell_dis_map[mcell]=dis;
	mcell_point_map[mcell]=p1;
      }
    }
  }
  
  return mcell_point_map;
}

double WCP::ToyPointCloud::get_closest_dis(WCP::Point& p){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,1);
  return sqrt(results.front().second);
}
std::pair<double, WCP::Point> WCP::ToyPointCloud::get_closest_point(WCP::Point& p){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,1);
  Point p1;
  int index = results.front().first;
  p1.x =  cloud.pts[index].x;
  p1.y =  cloud.pts[index].y;
  p1.z =  cloud.pts[index].z;
  return std::make_pair(sqrt(results.front().second),p1);
}



std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> WCP::ToyPointCloud::get_closest_points(WCP::Point& p, int N){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,N);
  std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> points;
  for (auto it = results.begin(); it!= results.end(); it++){
    size_t index = (*it).first;
    Point p;
    p.x = cloud.pts[index].x;
    p.y = cloud.pts[index].y;
    p.z = cloud.pts[index].z;
    SlimMergeGeomCell *mcell = cloud.pts[index].mcell;
    points.push_back(std::make_pair(mcell,p));
  }
  return points;
}
std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> WCP::ToyPointCloud::get_closest_points(WCP::Point& p, double search_radius){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,search_radius);
  std::vector<std::pair<WCP::SlimMergeGeomCell*,Point>> points;
  for (auto it = results.begin(); it!= results.end(); it++){
    size_t index = (*it).first;
    Point p;
    p.x = cloud.pts[index].x;
    p.y = cloud.pts[index].y;
    p.z = cloud.pts[index].z;
    SlimMergeGeomCell *mcell = cloud.pts[index].mcell;
    points.push_back(std::make_pair(mcell,p));
  }

  return points;
}



std::vector<WCPointCloud<double>::WCPoint> WCP::ToyPointCloud::get_hull(){
  quickhull::QuickHull<float> qh;
  std::vector<quickhull::Vector3<float>> pc;
  // std::cout << cloud.pts.size() << std::endl;
  for (size_t i=0;i!=cloud.pts.size();i++){
    // std::cout << cloud.pts.at(i).x << " " << cloud.pts.at(i).y << " " << cloud.pts.at(i).z << std::endl;
    pc.emplace_back(cloud.pts.at(i).x,cloud.pts.at(i).y,cloud.pts.at(i).z);
  }
  quickhull::ConvexHull<float> hull = qh.getConvexHull(pc,false,true);
  std::set<int> indices;
  
  for (size_t i=0;i!=hull.getIndexBuffer().size();i++){
    indices.insert(hull.getIndexBuffer().at(i));
  }
  
  std::vector<WCPointCloud<double>::WCPoint> results;
  //std::copy(indices.begin(),indices.end(),results.begin());
  for (auto it = indices.begin(); it!=indices.end(); it++){
    //std::cout << pc.at(*it).x << " " << pc.at(*it).y <<  " " << pc.at(*it).z << std::endl;
  // Point p;
  // p.x = cloud.pts.at(*it).x;
  // p.y = cloud.pts.at(*it).y;
  // p.z = cloud.pts.at(*it).z;
  // SlimMergeGeomCell *mcell = cloud.pts.at(*it).mcell;
  // results.push_back(std::make_pair(mcell,p));
  //results.push_back(*it);
    results.push_back(cloud.pts.at(*it));
  }
  return results;
}
