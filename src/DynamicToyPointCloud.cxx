#include "WireCellData/DynamicToyPointCloud.h"
#include "TH2F.h"

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



std::vector<std::pair<size_t,double>> WireCell::DynamicToyPointCloud::get_closest_index(WireCell::Point& p, int N){
  const size_t num_results = N;
  nanoflann::KNNResultSet<double> resultSet(num_results);
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  resultSet.init(&ret_index[0],&out_dist_sqr[0]);
  
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;  
  index->findNeighbors(resultSet,query_pt, nanoflann::SearchParams(10));
  
  std::vector<std::pair<size_t,double>> results(N);
  for(size_t i=0; i!=N; i++){
    results.at(i) = std::make_pair(ret_index.at(i),out_dist_sqr.at(i));
  }
  
  return results;
}



std::vector<std::pair<size_t,double>> WireCell::DynamicToyPointCloud::get_closest_2d_index(double x, double y, int N, int plane){
  const size_t num_results = N;
  nanoflann::KNNResultSet<double> resultSet(num_results);
  std::vector<size_t> ret_index(N);
  std::vector<double> out_dist_sqr(N);
  resultSet.init(&ret_index[0],&out_dist_sqr[0]);
  
  double query_pt[2];
  query_pt[0] = x;
  query_pt[1] = y;
  if (plane==0){
    index_u->findNeighbors(resultSet,query_pt, nanoflann::SearchParams(10));
  }else if (plane==1){
    index_v->findNeighbors(resultSet,query_pt, nanoflann::SearchParams(10));
  }else{
    index_w->findNeighbors(resultSet,query_pt, nanoflann::SearchParams(10));
  }
 
  
  std::vector<std::pair<size_t,double>> results(N);
  for(size_t i=0; i!=N; i++){
    results.at(i) = std::make_pair(ret_index.at(i),out_dist_sqr.at(i));
  }
  
  return results;
}


std::vector<std::pair<size_t,double>> WireCell::DynamicToyPointCloud::get_closest_index(WireCell::Point& p, double search_radius){
  double query_pt[3];
  query_pt[0] = p.x;
  query_pt[1] = p.y;
  query_pt[2] = p.z;
  std::vector<std::pair<size_t, double> > indices_dists;
  nanoflann::RadiusResultSet<double, size_t> resultSet(search_radius * search_radius, indices_dists);
  index->findNeighbors(resultSet, query_pt, nanoflann::SearchParams());
  
  
  return indices_dists;
}


std::vector<std::pair<size_t,double>> WireCell::DynamicToyPointCloud::get_closest_2d_index(double x, double y, double search_radius, int plane){
  double query_pt[2];
  query_pt[0] = x;
  query_pt[1] = y;
  std::vector<std::pair<size_t, double> > indices_dists;
  nanoflann::RadiusResultSet<double, size_t> resultSet(search_radius * search_radius, indices_dists);
  
  if (plane ==0){
    index_u->findNeighbors(resultSet, query_pt, nanoflann::SearchParams());
  }else if (plane==1){
    index_v->findNeighbors(resultSet, query_pt, nanoflann::SearchParams());
  }else{
    index_w->findNeighbors(resultSet, query_pt, nanoflann::SearchParams());
  }
  return indices_dists;
}


PR3DCluster* WireCell::DynamicToyPointCloud::get_cluster(int index){
  if (index <0 || index >= vec_index_cluster.size()){
    return 0;
  }else{
    return vec_index_cluster.at(index);
  }
}


void WireCell::DynamicToyPointCloud::AddPoints(PR3DCluster* cluster, Point& p_test, TVector3& dir, double range, double step, double angle){
  
  size_t current_size = cloud.pts.size();
  dir.SetMag(1);
  
  int num_points = int(range/(step))+1;
  double dis_seg = range/num_points;

  cloud.pts.resize(current_size + num_points);
  cloud_u.pts.resize(current_size + num_points);
  cloud_v.pts.resize(current_size + num_points);
  cloud_w.pts.resize(current_size + num_points);
  vec_index_cluster.resize(current_size + num_points);
  
  for (int k=0;k!=num_points;k++){
    double dis_cut = std::max(2.4*units::cm,k*dis_seg*sin(angle/180.*3.1415926));
    
    vec_index_cluster.at(current_size+k) = cluster;
    
    cloud.pts[current_size+k].x = p_test.x + k * dir.X() * dis_seg;
    cloud.pts[current_size+k].y = p_test.y + k * dir.Y() * dis_seg;
    cloud.pts[current_size+k].z = p_test.z + k * dir.Z() * dis_seg;
    cloud.pts[current_size+k].index_u = int(dis_cut); // hack ... 
    cloud.pts[current_size+k].index_v = int(dis_cut);
    cloud.pts[current_size+k].index_w = int(dis_cut); 
    cloud.pts[current_size+k].mcell = 0;
    cloud.pts[current_size+k].index = current_size+k;
      
    cloud_u.pts[current_size+k].x = cloud.pts[current_size+k].x;
    cloud_u.pts[current_size+k].y = cos(angle_u) * (cloud.pts[current_size+k].z) - sin(angle_u) * (cloud.pts[current_size+k].y);
    cloud_u.pts[current_size+k].index = current_size+k;
    
    cloud_v.pts[current_size+k].x = cloud.pts[current_size+k].x;
    cloud_v.pts[current_size+k].y = cos(angle_v) * (cloud.pts[current_size+k].z) - sin(angle_v) * (cloud.pts[current_size+k].y);
    cloud_v.pts[current_size+k].index = current_size+k;
    
    cloud_w.pts[current_size+k].x = cloud.pts[current_size+k].x;
    cloud_w.pts[current_size+k].y = cos(angle_w) * (cloud.pts[current_size+k].z) - sin(angle_w) * (cloud.pts[current_size+k].y);
    cloud_w.pts[current_size+k].index = current_size+k;
  }

  index->addPoints(current_size, current_size+num_points-1);
  index_u->addPoints(current_size, current_size+num_points-1);
  index_v->addPoints(current_size, current_size+num_points-1);
  index_w->addPoints(current_size, current_size+num_points-1);
  
}

void WireCell::DynamicToyPointCloud::AddPoints(PR3DCluster* cluster, int flag, double step){
  size_t current_size = cloud.pts.size();
  
  if (flag==0){
    // add actual points in
    WireCell::WCPointCloud<double>& pcloud = cluster->get_point_cloud()->get_cloud();
    WireCell::WC2DPointCloud<double>& pcloud_u = cluster->get_point_cloud()->get_cloud_u();
    WireCell::WC2DPointCloud<double>& pcloud_v = cluster->get_point_cloud()->get_cloud_v();
    WireCell::WC2DPointCloud<double>& pcloud_w = cluster->get_point_cloud()->get_cloud_w();

    cloud.pts.resize(current_size + pcloud.pts.size());
    cloud_u.pts.resize(current_size + pcloud.pts.size());
    cloud_v.pts.resize(current_size + pcloud.pts.size());
    cloud_w.pts.resize(current_size + pcloud.pts.size());
    vec_index_cluster.resize(current_size + pcloud.pts.size());
    
    for (size_t i=0;i!=pcloud.pts.size();i++){
      vec_index_cluster.at(current_size+i) = cluster;

      cloud.pts[current_size+i].x = pcloud.pts.at(i).x;
      cloud.pts[current_size+i].y = pcloud.pts.at(i).y;
      cloud.pts[current_size+i].z = pcloud.pts.at(i).z;
      cloud.pts[current_size+i].index_u = pcloud.pts.at(i).index_u;
      cloud.pts[current_size+i].index_v = pcloud.pts.at(i).index_v;
      cloud.pts[current_size+i].index_w = pcloud.pts.at(i).index_w;
      cloud.pts[current_size+i].mcell = pcloud.pts.at(i).mcell;
      cloud.pts[current_size+i].index = current_size+i;

      cloud_u.pts[current_size+i].x = pcloud_u.pts.at(i).x;
      cloud_u.pts[current_size+i].y = pcloud_u.pts.at(i).y;
      cloud_u.pts[current_size+i].index = current_size+i;
      
      cloud_v.pts[current_size+i].x = pcloud_v.pts.at(i).x;
      cloud_v.pts[current_size+i].y = pcloud_v.pts.at(i).y;
      cloud_v.pts[current_size+i].index = current_size+i;
      
      cloud_w.pts[current_size+i].x = pcloud_w.pts.at(i).x;
      cloud_w.pts[current_size+i].y = pcloud_w.pts.at(i).y;
      cloud_w.pts[current_size+i].index = current_size+i;
    }
    if (pcloud.pts.size()>0){
      index->addPoints(current_size, current_size+pcloud.pts.size()-1);
      index_u->addPoints(current_size, current_size+pcloud.pts.size()-1);
      index_v->addPoints(current_size, current_size+pcloud.pts.size()-1);
      index_w->addPoints(current_size, current_size+pcloud.pts.size()-1);
    }
  }else{
    // add skeleton points in
    std::list<WCPointCloud<double>::WCPoint>& path_wcps = cluster->get_path_wcps();

    PointVector pts;
    WCPointCloud<double>::WCPoint prev_wcp = path_wcps.front();
    for (auto it = path_wcps.begin(); it!=path_wcps.end();it++){
      
      double dis = sqrt(pow((*it).x - prev_wcp.x,2) + pow((*it).y - prev_wcp.y,2) + pow((*it).z - prev_wcp.z,2));
      if (dis <=step){
	Point current_pt((*it).x,(*it).y,(*it).z);
	//	std::cout << current_pt.x/units::cm << " " << current_pt.y/units::cm << " " << current_pt.z/units::cm << std::endl;
	pts.push_back(current_pt);
      }else{
	int num_points = int(dis/(step))+1;
	double dis_seg = dis/num_points;
	for (int k=0;k!=num_points;k++){
	  Point current_pt(prev_wcp.x + (k+1.)/num_points*((*it).x - prev_wcp.x),
			   prev_wcp.y + (k+1.)/num_points*((*it).y - prev_wcp.y),
			   prev_wcp.z + (k+1.)/num_points*((*it).z - prev_wcp.z));
	  //  std::cout << current_pt.x/units::cm << " " << current_pt.y/units::cm << " " << current_pt.z/units::cm << std::endl;

	  pts.push_back(current_pt);
	  
	}
      }
      prev_wcp = (*it);
    }

    
    cloud.pts.resize(current_size + pts.size());
    cloud_u.pts.resize(current_size + pts.size());
    cloud_v.pts.resize(current_size + pts.size());
    cloud_w.pts.resize(current_size + pts.size());
    vec_index_cluster.resize(current_size + pts.size());
    int i = 0;
    for (auto it = pts.begin(); it!=pts.end();it++){
      vec_index_cluster.at(current_size+i) = cluster;

      cloud.pts[current_size+i].x = (*it).x;
      cloud.pts[current_size+i].y = (*it).y;
      cloud.pts[current_size+i].z = (*it).z;
      cloud.pts[current_size+i].index_u = 2.4*units::cm;
      cloud.pts[current_size+i].index_v = 2.4*units::cm;
      cloud.pts[current_size+i].index_w = 2.4*units::cm;
      cloud.pts[current_size+i].mcell = 0;
      cloud.pts[current_size+i].index = current_size+i;

      cloud_u.pts[current_size+i].x = (*it).x;
      cloud_u.pts[current_size+i].y = cos(angle_u) * (*it).z - sin(angle_u) * (*it).y;
      cloud_u.pts[current_size+i].index = current_size+i;
      
      cloud_v.pts[current_size+i].x = (*it).x;
      cloud_v.pts[current_size+i].y = cos(angle_v) * (*it).z - sin(angle_v) * (*it).y;
      cloud_v.pts[current_size+i].index = current_size+i;
      
      cloud_w.pts[current_size+i].x = (*it).x;
      cloud_w.pts[current_size+i].y = cos(angle_w) * (*it).z - sin(angle_w) * (*it).y;
      cloud_w.pts[current_size+i].index = current_size+i;
      
      i ++;
    }
    if (pts.size()>0){
      index->addPoints(current_size, current_size+pts.size()-1);
      index_u->addPoints(current_size, current_size+pts.size()-1);
      index_v->addPoints(current_size, current_size+pts.size()-1);
      index_w->addPoints(current_size, current_size+pts.size()-1);
    }
  }
}



std::tuple<double, PR3DCluster*, size_t>  WireCell::DynamicToyPointCloud::get_closest_point_info(WireCell::Point& p){
  std::vector<std::pair<size_t,double>> results = get_closest_index(p,1);
  return std::make_tuple(sqrt(results.front().second), vec_index_cluster.at(results.front().first), results.front().first);
}

std::vector<std::tuple<double, PR3DCluster*, size_t>> WireCell::DynamicToyPointCloud::get_2d_points_info(WireCell::Point& p, double radius, int plane){
  std::vector<std::pair<size_t,double>> results;
  double x,y;
  if (plane==0){
    x = p.x;
    y = cos(angle_u) * p.z - sin(angle_u) * p.y;
    results = get_closest_2d_index(x,y,radius,0);
  }else if (plane==1){
    x = p.x;
    y = cos(angle_v) * p.z - sin(angle_v) * p.y;
    results = get_closest_2d_index(x,y,radius,1);
  }else if (plane==2){
    x = p.x;
    y = cos(angle_w) * p.z - sin(angle_w) * p.y;
    results = get_closest_2d_index(x,y,radius,2);
  }
  std::vector<std::tuple<double, PR3DCluster*, size_t>> return_results;

  for (size_t i=0;i!=results.size();i++){
    return_results.push_back(std::make_tuple(sqrt(results.at(i).second), vec_index_cluster.at(results.at(i).first), results.at(i).first));
  }
  
  
  return return_results;
}

std::tuple<double, PR3DCluster*, size_t> WireCell::DynamicToyPointCloud::get_closest_2d_point_info(WireCell::Point& p, int plane){
  std::vector<std::pair<size_t,double>> results;
  double x,y;
  if (plane==0){
    x = p.x;
    y = cos(angle_u) * p.z - sin(angle_u) * p.y;
    results = get_closest_2d_index(x,y,1,0);
  }else if (plane==1){
    x = p.x;
    y = cos(angle_v) * p.z - sin(angle_v) * p.y;
    results = get_closest_2d_index(x,y,1,1);
  }else if (plane==2){
    x = p.x;
    y = cos(angle_w) * p.z - sin(angle_w) * p.y;
    results = get_closest_2d_index(x,y,1,2);
  }
  return std::make_tuple(sqrt(results.front().second), vec_index_cluster.at(results.front().first), results.front().first);
}

std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> WireCell::DynamicToyPointCloud::get_closest_points(WireCell::Point& p, double search_radius){
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




std::pair<double,double> WireCell::DynamicToyPointCloud::HoughTrans(Point&p , double dis){
  double theta, phi;
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  
  std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> pts = get_closest_points(p,dis);

  // std::cout << "Num " <<  pts.size() << std::endl;
    
  double x,y,z,q;
  for (size_t i=0; i!=pts.size(); i++){
    x = pts.at(i).second.x;
    y = pts.at(i).second.y;
    z = pts.at(i).second.z;
    q = pts.at(i).first->get_q()/pts.at(i).first->get_sampling_points().size();
    if (q<=0) continue;
    double r = sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));

    //  for (int i1=0; i1!=5; i1++){
    //  for (int j1=0; j1!=5; j1++){
    //	for (int k1=0; k1!=5; k1++){
    TVector3 vec(x-x0 ,y-y0 ,z-z0 );
    if ( r < 10*units::cm){
      hough->Fill(vec.Theta(),vec.Phi(), q );
    }else{
      hough->Fill(vec.Theta(),vec.Phi(), q *pow(10*units::cm/r,2) );
    }
  }
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta = hough->GetXaxis()->GetBinCenter(a);
  phi = hough->GetYaxis()->GetBinCenter(b);

  // std::cout << hough->GetSum() << " " << hough->GetBinContent(a,b)<< std::endl;
  
  delete hough;
  return std::make_pair(theta,phi);
}

TVector3 WireCell::DynamicToyPointCloud::VHoughTrans(Point& p, double dis){
  double theta, phi;
  std::pair<double,double> angles_1 = HoughTrans(p,dis);
  theta = angles_1.first;
  phi = angles_1.second;
  TVector3 temp(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  return temp;
}

