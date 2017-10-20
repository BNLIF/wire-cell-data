#include "WireCellData/ToyPointCloud.h"

using namespace WireCell;

WireCell::ToyPointCloud::ToyPointCloud(){
  index = 0;
}

WireCell::ToyPointCloud::~ToyPointCloud(){
  if (index!=0)
    delete index;

  cloud.pts.clear();
}

void WireCell::ToyPointCloud::AddPoint(Point& p, SlimMergeGeomCell *mcell){
  WCPointCloud<double>::WCPoint point;
  point.x = p.x;
  point.y = p.y;
  point.z = p.z;
  point.mcell = mcell;
  cloud.pts.push_back(point);
}

void WireCell::ToyPointCloud::AddPoints(PointVector& ps, SlimMergeGeomCell *mcell){
  size_t current_size = cloud.pts.size();
  cloud.pts.resize(current_size + ps.size());
  for (size_t i=0;i!=ps.size();i++){
    cloud.pts[current_size+i].x = ps.at(i).x;
    cloud.pts[current_size+i].y = ps.at(i).y;
    cloud.pts[current_size+i].z = ps.at(i).z;
    cloud.pts[current_size+i].mcell = mcell;
  }
}

void WireCell::ToyPointCloud::build_kdtree_index(){
  if (index!=(my_kd_tree_t*)0){
    delete index;
  }
  index = new my_kd_tree_t(3 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  index->buildIndex();
}

std::vector<std::pair<size_t,double>> WireCell::ToyPointCloud::get_closest_index(WireCell::Point& p, int N){
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

std::vector<std::pair<size_t,double>> WireCell::ToyPointCloud::get_closest_index(WireCell::Point& p, double search_radius){
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;
  std::vector<std::pair<size_t,double> >   ret_matches;
  nanoflann::SearchParams params;
  const size_t nMatches = index->radiusSearch(&query_pt[0], search_radius*search_radius, ret_matches, params);
  
  return ret_matches;
}

std::map<WireCell::SlimMergeGeomCell*, Point> WireCell::ToyPointCloud::get_closest_mcell(WireCell::Point& p, int N){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,N);
  
  std::map<WireCell::SlimMergeGeomCell*, Point> mcell_point_map;
  std::map<WireCell::SlimMergeGeomCell*, double> mcell_dis_map;
  
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


std::map<WireCell::SlimMergeGeomCell*, Point> WireCell::ToyPointCloud::get_closest_mcell(WireCell::Point& p, double search_radius){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,search_radius);
  
  std::map<WireCell::SlimMergeGeomCell*, Point> mcell_point_map;
  std::map<WireCell::SlimMergeGeomCell*, double> mcell_dis_map;
  
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

std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> WireCell::ToyPointCloud::get_closest_points(WireCell::Point& p, int N){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,N);
  std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> points;
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
std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> WireCell::ToyPointCloud::get_closest_points(WireCell::Point& p, double search_radius){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,search_radius);
  std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> points;
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



std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> WireCell::ToyPointCloud::get_hull(){
  quickhull::QuickHull<float> qh;
  std::vector<quickhull::Vector3<float>> pc;
  for (size_t i=0;i!=cloud.pts.size();i++){
    pc.emplace_back(cloud.pts.at(i).x,cloud.pts.at(i).y,cloud.pts.at(i).z);
  }
  quickhull::ConvexHull<float> hull = qh.getConvexHull(pc,false,false);
  std::set<int> indices;
  
  for (size_t i=0;i!=hull.getIndexBuffer().size();i++){
    indices.insert(hull.getIndexBuffer().at(i));
  }
  
  std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> results;
  for (auto it = indices.begin(); it!=indices.end(); it++){
    Point p;
    p.x = cloud.pts.at(*it).x;
    p.y = cloud.pts.at(*it).y;
    p.z = cloud.pts.at(*it).z;
    SlimMergeGeomCell *mcell = cloud.pts.at(*it).mcell;
    results.push_back(std::make_pair(mcell,p));
  }
  return results;
}