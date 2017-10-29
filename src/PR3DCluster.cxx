#include "WireCellData/PR3DCluster.h"

#include "TMatrixDEigen.h"
#include "TH2F.h"

#include <boost/graph/connected_components.hpp>



using namespace WireCell;

PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
  point_cloud = 0;
  graph = 0;
}

PR3DCluster::~PR3DCluster(){
  if (point_cloud!=(ToyPointCloud*)0)
    delete point_cloud;
  if (graph!=(MCUGraph*)0)
    delete graph;
}

// void AddCell(SlimMergeGeomCell* mcell, int *time_slices, int ntime_slice){
// }

void PR3DCluster::Create_point_cloud(){
  if (point_cloud!=(ToyPointCloud*)0)
    return;
  
  point_cloud = new ToyPointCloud();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    PointVector pts = mcell->get_sampling_points();
    
    point_cloud->AddPoints(pts,mcell->get_sampling_points_wires(),mcell);
  }
  point_cloud->build_kdtree_index();
  //  std::cout << point_cloud->get_num_points() << std::endl;

}


void PR3DCluster::Create_graph(){
  if (point_cloud==(ToyPointCloud*)0)
    Create_point_cloud();

  // create Graph ...
  const int N = point_cloud->get_num_points();
  graph = new MCUGraph(N);

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_uindex_wcps;
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_vindex_wcps;
  std::map<SlimMergeGeomCell*, std::map<int, std::set<int>>> map_mcell_windex_wcps;
  
  //point_cloud->Print();
  
  //create all vertices
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::map<int, std::set<int>> map_uindex_wcps;
    std::map<int, std::set<int>> map_vindex_wcps;
    std::map<int, std::set<int>> map_windex_wcps;
    std::vector<int>& wcps = point_cloud->get_mcell_indices(mcell);
    for (auto it1 = wcps.begin(); it1!=wcps.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp = cloud.pts[*it1]; 
      // int index = wcp.index;
      // std::cout << index << " " << wcp.x << " " << wcp.y << " " << wcp.z << " " << wcp.index_u << " " << wcp.index_v << " " << wcp.index_w << std::endl;
      
      auto v = vertex(wcp.index, *graph); // retrieve vertex descriptor
      (*graph)[v].index = wcp.index;
      if (map_uindex_wcps.find(wcp.index_u)==map_uindex_wcps.end()){
   	std::set<int> wcps;
   	wcps.insert(wcp.index);
   	map_uindex_wcps[wcp.index_u] = wcps;
      }else{
   	map_uindex_wcps[wcp.index_u].insert(wcp.index);
      }
      
      if (map_vindex_wcps.find(wcp.index_v)==map_vindex_wcps.end()){
  	std::set<int> wcps;
  	wcps.insert(wcp.index);
  	map_vindex_wcps[wcp.index_v] = wcps;
      }else{
  	map_vindex_wcps[wcp.index_v].insert(wcp.index);
      }

      if (map_windex_wcps.find(wcp.index_w)==map_windex_wcps.end()){
  	std::set<int> wcps;
  	wcps.insert(wcp.index);
  	map_windex_wcps[wcp.index_w] = wcps;
      }else{
  	map_windex_wcps[wcp.index_w].insert(wcp.index);
      }
      
      
    }
    map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
    map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
    map_mcell_windex_wcps[mcell] = map_windex_wcps;
  }
  
  // create graph for points inside the same mcell
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    std::vector<int>& wcps = point_cloud->get_mcell_indices(mcell);
    int max_wire_interval = mcell->get_max_wire_interval();
    int min_wire_interval = mcell->get_min_wire_interval();
    std::map<int, std::set<int>>* map_max_index_wcps;
    std::map<int, std::set<int>>* map_min_index_wcps;
    if (mcell->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
    }else if (mcell->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell];
    }
    if (mcell->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
    }else if (mcell->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell];
    }
    
    for (auto it1 = wcps.begin(); it1!=wcps.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }

      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }
      // std::cout << max_wcps_set.size() << " " << min_wcps_set.size() << std::endl;
      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
   	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
   	  std::set<int> common_set;
   	  set_intersection((*it2)->begin(), (*it2)->end(), (*it3)->begin(), (*it3)->end(),std::inserter(common_set,common_set.begin()));
  	  for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
   	    WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
   	    if (wcp2.index != wcp1.index){
   	      int index2 = wcp2.index;
   	      //  std::cout << index1 << " " << index2 << std::endl;
   	      // add edge ...
   	      auto edge = add_edge(index1,index2,*graph);
   	      if (edge.second){
   		(*graph)[edge.first].dist = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
		//	std::cout << wcp1.x << " " << wcp1.y << " " << wcp1.z << " " << wcp1.index_u << " " << wcp1.index_v << " " << wcp1.index_w << " " << wcp2.index_u << " " << wcp2.index_v << " " << wcp2.index_w << std::endl;
   	      }
   	    }
   	  }
   	}
      }
    }
  }

  
  
  
  std::vector<int> time_slices;
  for (auto it1 = time_cells_set_map.begin(); it1!=time_cells_set_map.end(); it1++){
    time_slices.push_back((*it1).first);
  }

  std::vector<std::pair<SlimMergeGeomCell*,SlimMergeGeomCell*>> connected_mcells;
  
  for (size_t i=0; i!= time_slices.size(); i++){
    SMGCSet& mcells_set = time_cells_set_map[time_slices.at(i)];
    
    // create graph for points in mcell inside the same time slice
    if (mcells_set.size()>=2){
      for (auto it2 = mcells_set.begin(); it2!=mcells_set.end();it2++){
  	SlimMergeGeomCell *mcell1 = *it2;
  	auto it2p = it2;
  	if (it2p!=mcells_set.end()){
  	  it2p++;
  	  for (auto it3 = it2p; it3!=mcells_set.end(); it3++){
  	    SlimMergeGeomCell *mcell2 = *(it3);
  	    //std::cout << mcell1 << " " << mcell2 << " " << mcell1->Overlap_fast(mcell2,2) << std::endl;
  	    connected_mcells.push_back(std::make_pair(mcell1,mcell2));
  	  }
  	}
      }
    }
    // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
    std::vector<SMGCSet> vec_mcells_set;
    if (i+1 < time_slices.size()){
      if (time_slices.at(i+1)-time_slices.at(i)==1){
	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
	if (i+2 < time_slices.size())
	  if (time_slices.at(i+2)-time_slices.at(i)==2)
	    vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+2)]);
      }else if (time_slices.at(i+1) - time_slices.at(i)==2){
	vec_mcells_set.push_back(time_cells_set_map[time_slices.at(i+1)]);
      }
    }
    bool flag = false;
    for (size_t j=0; j!=vec_mcells_set.size(); j++){
      if (flag) break;
      SMGCSet& next_mcells_set = vec_mcells_set.at(j);
      for (auto it1 = mcells_set.begin(); it1!= mcells_set.end(); it1++){
	SlimMergeGeomCell *mcell1 = (*it1);
	for (auto it2 = next_mcells_set.begin(); it2!=next_mcells_set.end(); it2++){
	  SlimMergeGeomCell *mcell2 = (*it2);
	  if (mcell1->Overlap_fast(mcell2,2)){
	    flag = true;
	    connected_mcells.push_back(std::make_pair(mcell1,mcell2));
	  }
	}
      }
    }
  }
  
  // establish edge ... 
  // std::cout << connected_mcells.size() << std::endl;
  for (auto it = connected_mcells.begin(); it!= connected_mcells.end(); it++){
    SlimMergeGeomCell *mcell1 = (*it).first;
    SlimMergeGeomCell *mcell2 = (*it).second;

    std::vector<int>& wcps1 = point_cloud->get_mcell_indices(mcell1);
    std::vector<int>& wcps2 = point_cloud->get_mcell_indices(mcell2);

    // test 2 against 1 ... 
    int max_wire_interval = mcell1->get_max_wire_interval();
    int min_wire_interval = mcell1->get_min_wire_interval();
    std::map<int, std::set<int>>* map_max_index_wcps;
    std::map<int, std::set<int>>* map_min_index_wcps;
    
    if (mcell1->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell2];
    }else if (mcell1->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell2];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell2];
    }
    if (mcell1->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell2];
    }else if (mcell1->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell2];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell2];
    }

    for (auto it1 = wcps1.begin(); it1!=wcps1.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell1->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell1->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell1->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell1->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }
      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }

      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
   	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
   	  std::set<int> common_set;
   	  set_intersection((*it2)->begin(), (*it2)->end(), (*it3)->begin(), (*it3)->end(),std::inserter(common_set,common_set.begin()));
  	  for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
   	    WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
   	    if (wcp2.index != wcp1.index){
   	      int index2 = wcp2.index;
	      auto edge = add_edge(index1,index2,*graph);
   	      if (edge.second){
   		(*graph)[edge.first].dist = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	      }
   	    }
   	  }
   	}
      }
      
    }


    // test 1 against 2 ...
    max_wire_interval = mcell2->get_max_wire_interval();
    min_wire_interval = mcell2->get_min_wire_interval();
    if (mcell2->get_max_wire_type()==WirePlaneType_t(0)){
      map_max_index_wcps = &map_mcell_uindex_wcps[mcell1];
    }else if (mcell2->get_max_wire_type()==WirePlaneType_t(1)){
      map_max_index_wcps = &map_mcell_vindex_wcps[mcell1];
    }else{
      map_max_index_wcps = &map_mcell_windex_wcps[mcell1];
    }
    if (mcell2->get_min_wire_type()==WirePlaneType_t(0)){
      map_min_index_wcps = &map_mcell_uindex_wcps[mcell1];
    }else if (mcell2->get_min_wire_type()==WirePlaneType_t(1)){
      map_min_index_wcps = &map_mcell_vindex_wcps[mcell1];
    }else{
      map_min_index_wcps = &map_mcell_windex_wcps[mcell1];
    }
    for (auto it1 = wcps2.begin(); it1!=wcps2.end(); it1++){
      WCPointCloud<double>::WCPoint& wcp1 = cloud.pts[*it1];
      int index1 = wcp1.index;
      int index_max_wire;
      int index_min_wire;
      if (mcell2->get_max_wire_type()==WirePlaneType_t(0)){
  	index_max_wire = wcp1.index_u;
      }else if (mcell2->get_max_wire_type()==WirePlaneType_t(1)){
  	index_max_wire = wcp1.index_v;
      }else{
  	index_max_wire = wcp1.index_w;
      }
      if (mcell2->get_min_wire_type()==WirePlaneType_t(0)){
  	index_min_wire = wcp1.index_u;
      }else if (mcell2->get_min_wire_type()==WirePlaneType_t(1)){
  	index_min_wire = wcp1.index_v;
      }else{
  	index_min_wire = wcp1.index_w;
      }
      std::vector<std::set<int>*> max_wcps_set;
      std::vector<std::set<int>*> min_wcps_set;
      // go through the first map and find the ones satisfying the condition
      for (auto it2 = map_max_index_wcps->begin(); it2!=map_max_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_max_wire)<=max_wire_interval){
   	  max_wcps_set.push_back(&(it2->second));
   	}
      }
      // go through the second map and find the ones satisfying the condition
      for (auto it2 = map_min_index_wcps->begin(); it2!=map_min_index_wcps->end(); it2++){
   	if (fabs(it2->first - index_min_wire)<=min_wire_interval){
   	  min_wcps_set.push_back(&(it2->second));
   	}
      }

      for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
   	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
   	  std::set<int> common_set;
   	  set_intersection((*it2)->begin(), (*it2)->end(), (*it3)->begin(), (*it3)->end(),std::inserter(common_set,common_set.begin()));
  	  for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
   	    WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
   	    if (wcp2.index != wcp1.index){
   	      int index2 = wcp2.index;
	      auto edge = add_edge(index1,index2,*graph);
   	      if (edge.second){
   		(*graph)[edge.first].dist = sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2));
	      }
   	    }
   	  }
   	}
      }
    }
  }
  
  // get the connected components from the graph
  {
    std::vector<int> component(num_vertices(*graph));
    const int num = connected_components(*graph,&component[0]);

    if (num >1){
      //for separated kd tree to find the closest points between disconnected components,
      std::vector<ToyPointCloud*> pt_clouds;
      for (int j=0;j!=num;j++){
	ToyPointCloud *pt_cloud = new ToyPointCloud();
	pt_clouds.push_back(pt_cloud);
      }
      
      std::vector<int>::size_type i;
      for (i=0;i!=component.size(); ++i){
	pt_clouds.at(component[i])->AddPoint(cloud.pts[i]);
	//   std::cout << "Vertex " << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " << cloud.pts[i].mcell << " " << cloud.pts[i].mcell->GetTimeSlice()  << " is in component " << component[i] << std::endl;
      }
      for (int j=0;j!=num;j++){
	pt_clouds.at(j)->build_kdtree_index();
      }
      // connect these graphs according to closest distance some how ...
      
      
      for (int i=0;i!=num;i++){
	delete pt_clouds.at(i);
      }
    }

    
  }

  



  
}

Point PR3DCluster::calc_ave_pos(Point& p, double dis){
  std::map<WireCell::SlimMergeGeomCell*, Point> pts = point_cloud->get_closest_mcell(p,dis);
  Point pt(0,0,0);
  double charge = 0;
  //std::cout << pts.size() << std::endl;
  for (auto it = pts.begin(); it!= pts.end(); it++){
    SlimMergeGeomCell *mcell = (*it).first;
    Point pc = mcell->center();
    double q = mcell->get_q();
    pt.x += pc.x * q;
    pt.y += pc.y * q;
    pt.z += pc.z * q;
    charge += q;
    // std::cout << pc.x/units::cm << " " << pc.y/units::cm << " " << pc.z/units::cm << " " << q << " " <<
    //  sqrt(pow(pc.x-p.x,2)+pow(pc.y-p.y,2)+pow(pc.z-p.z,2))/units::cm << std::endl;
  }
  if (charge!=0){
    pt.x/=charge;
    pt.y/=charge;
    pt.z/=charge;
  }
  return pt;
    
}

std::pair<SlimMergeGeomCell*, Point> PR3DCluster::get_closest_point_mcell(Point& p_test){
  std::map<SlimMergeGeomCell*,Point> pts = point_cloud->get_closest_mcell(p_test,1);
  return std::make_pair((*pts.begin()).first,(*pts.begin()).second);
}

TVector3 PR3DCluster::calc_PCA_dir(Point& p, double dis){
  std::map<WireCell::SlimMergeGeomCell*, Point> pts = point_cloud->get_closest_mcell(p,dis);
  Point center1(0,0,0);
  // double charge=0;
  // for (auto it = pts.begin(); it!= pts.end(); it++){
  //   SlimMergeGeomCell *mcell = (*it).first;
  //   Point point = mcell->center();
  //   double q = mcell->get_q();
  //   center1.x += point.x * q;
  //   center1.y += point.y * q;
  //   center1.z += point.z * q;
  //   charge+=q;
  // }
  // center1.x/=charge;
  // center1.y/=charge;
  // center1.z/=charge;
  center1.x = p.x;
  center1.y = p.y;
  center1.z = p.z;

  TMatrixD cov_matrix(3,3);
  for (int i=0;i!=3;i++){
    for (int j=i;j!=3;j++){
      cov_matrix(i,j)=0;
      for (auto it = pts.begin(); it!= pts.end(); it++){
	SlimMergeGeomCell *mcell = (*it).first;
	double q = mcell->get_q();
	PointVector ps = mcell->get_sampling_points();
	for (int k=0;k!=ps.size();k++){
	  if (i==0 && j==0){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).x - center.x);//*q*q/ps.size()/ps.size();
	  }else if (i==0 && j==1){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).y - center.y);//*q*q/ps.size()/ps.size();
	  }else if (i==0 && j==2){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).z - center.z);//*q*q/ps.size()/ps.size();
	  }else if (i==1 && j==1){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).y - center.y);//*q*q/ps.size()/ps.size();
	  }else if (i==1 && j==2){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).z - center.z);//*q*q/ps.size()/ps.size();
	  }else if (i==2 && j==2){
	    cov_matrix(i,j) += (ps.at(k).z - center.z) * (ps.at(k).z - center.z);//*q*q/ps.size()/ps.size();
	  }
	}
      }
    }
  }

  cov_matrix(1,0) = cov_matrix(0,1);
  cov_matrix(2,0) = cov_matrix(0,2);
  cov_matrix(2,1) = cov_matrix(1,2);
  
  TMatrixDEigen eigen(cov_matrix);
  TMatrixD eigen_values = eigen.GetEigenValues();
  TMatrixD eigen_vectors = eigen.GetEigenVectors();
  TVector3 dir(eigen_vectors(0,0)/sqrt(eigen_vectors(0,0)*eigen_vectors(0,0) + eigen_vectors(1,0)*eigen_vectors(1,0) + eigen_vectors(2,0)*eigen_vectors(2,0)),
	       eigen_vectors(1,0)/sqrt(eigen_vectors(0,0)*eigen_vectors(0,0) + eigen_vectors(1,0)*eigen_vectors(1,0) + eigen_vectors(2,0)*eigen_vectors(2,0)),
	       eigen_vectors(2,0)/sqrt(eigen_vectors(0,0)*eigen_vectors(0,0) + eigen_vectors(1,0)*eigen_vectors(1,0) + eigen_vectors(2,0)*eigen_vectors(2,0)));
  return dir;
  
}

TVector3 PR3DCluster::calc_dir(Point& p_test, Point& p, double dis){
  std::map<WireCell::SlimMergeGeomCell*, Point> pts = point_cloud->get_closest_mcell(p,dis);
  TVector3 dir(0,0,0);
  for (auto it = pts.begin(); it!= pts.end(); it++){
    SlimMergeGeomCell *mcell = (*it).first;
    Point point = mcell->center();
    double q = mcell->get_q();
    TVector3 dir1(point.x-p_test.x,point.y-p_test.y,point.z-p_test.z);
    dir += dir1 * q;
  }
  if (dir.Mag()!=0)
    dir.SetMag(1);
  return dir;
}

TVector3 PR3DCluster::VHoughTrans(Point&p, double dis){
  double theta, phi;
  std::pair<double,double> angles_1 = HoughTrans(p,dis);
  theta = angles_1.first;
  phi = angles_1.second;
  TVector3 temp(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
  return temp;
}

std::pair<double,double> PR3DCluster::HoughTrans(Point&p , double dis){
  double theta, phi;
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  
  std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>>pts = point_cloud->get_closest_points(p,dis);

  // std::cout << "Num " <<  pts.size() << std::endl;
    
  double x,y,z,q;
  for (size_t i=0; i!=pts.size(); i++){
    x = pts.at(i).second.x;
    y = pts.at(i).second.y;
    z = pts.at(i).second.z;
    q = pts.at(i).first->get_q()/pts.at(i).first->get_sampling_points().size();
    if (q<=0) continue;

    //  for (int i1=0; i1!=5; i1++){
    //  for (int j1=0; j1!=5; j1++){
    //	for (int k1=0; k1!=5; k1++){
    TVector3 vec(x-x0 ,y-y0 ,z-z0 );
    hough->Fill(vec.Theta(),vec.Phi(), q );
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


void PR3DCluster::Calc_PCA(){
  center.x=0; center.y=0; center.z=0;
  int nsum = 0;
  for (auto it = mcells.begin(); it!=mcells.end();it++){
    PointVector ps = (*it)->get_sampling_points();
    for (int k=0;k!=ps.size();k++){
      center.x += ps.at(k).x;
      center.y += ps.at(k).y;
      center.z += ps.at(k).z;
      nsum ++;
    }
  }

  for (int i=0;i!=3;i++){
    PCA_axis[i].x = 0;
    PCA_axis[i].y = 0;
    PCA_axis[i].z = 0;
  }
  
  if (nsum>=3){
    center.x /=nsum;
    center.y /=nsum;
    center.z /=nsum;
  }else{
    return;
  }
  TMatrixD cov_matrix(3,3);

  for (int i=0;i!=3;i++){
    for (int j=i;j!=3;j++){
      cov_matrix(i,j)=0;
      for (auto it = mcells.begin(); it!=mcells.end();it++){
	PointVector ps = (*it)->get_sampling_points();
	for (int k=0;k!=ps.size();k++){
	  if (i==0 && j==0){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).x - center.x);
	  }else if (i==0 && j==1){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).y - center.y);
	  }else if (i==0 && j==2){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).z - center.z);
	  }else if (i==1 && j==1){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).y - center.y);
	  }else if (i==1 && j==2){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).z - center.z);
	  }else if (i==2 && j==2){
	    cov_matrix(i,j) += (ps.at(k).z - center.z) * (ps.at(k).z - center.z);
	  }
	}
      }
    }
  }
  cov_matrix(1,0) = cov_matrix(0,1);
  cov_matrix(2,0) = cov_matrix(0,2);
  cov_matrix(2,1) = cov_matrix(1,2);
  
  TMatrixDEigen eigen(cov_matrix);
  TMatrixD eigen_values = eigen.GetEigenValues();
  TMatrixD eigen_vectors = eigen.GetEigenVectors();

  // std::cout << eigen_values(0,0) << " " << eigen_values(1,1) << " " << eigen_values(2,2) << std::endl;
  //std::cout << eigen_vectors(0,0) << " " << eigen_vectors(0,1) << " " << eigen_vectors(0,2) << std::endl;
  for (int i=0;i!=3;i++){
    PCA_axis[i].x = eigen_vectors(0,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].y = eigen_vectors(1,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].z = eigen_vectors(2,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));

  }
  
  
  //std::cout << mean_x << " " << mean_y << " " << mean_z << std::endl;
}


void PR3DCluster::AddCell(SlimMergeGeomCell* mcell, int time_slice){
  if (cell_times_set_map.find(mcell)==cell_times_set_map.end()){
    std::set<int> times;
    times.insert(time_slice);
    cell_times_set_map[mcell]=times;
    mcells.push_back(mcell);
  }else{
    std::set<int>& times = cell_times_set_map[mcell];
    //if (find(times.begin(),times.end(),time_slice)==times.end()){
    times.insert(time_slice);
    // }
  }
  
  if (time_cells_set_map.find(time_slice)==time_cells_set_map.end()){
    SMGCSet mcells_1;
    mcells_1.insert(mcell);
    time_cells_set_map[time_slice] = mcells_1;
  }else{
    SMGCSet& mcells_1 = time_cells_set_map[time_slice];
    //if (find(mcells_1.begin(),mcells_1.end(), mcell) == mcells_1.end()){
    mcells_1.insert(mcell);
    //}
  }
}

void PR3DCluster::Remove_duplicated_mcells(){
  SMGCSelection to_be_removed;
  SMGCSelection to_be_saved;
  std::map<SlimMergeGeomCell*,SlimMergeGeomCell*> mcell_map;
  for (auto it = mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    bool save = true;
    for (auto it1=to_be_saved.begin(); it1!=to_be_saved.end(); it1++){
      SlimMergeGeomCell *mcell1 = *it1;
      if (mcell->IsSame(mcell1)){
	to_be_removed.push_back(mcell);
	mcell_map[mcell] = mcell1;
	save = false;
	break;
      }
    }
    if(save) to_be_saved.push_back(mcell);
  }
  // Now merge ...
  for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
    SlimMergeGeomCell *mcell = *it; // to be removed
    SlimMergeGeomCell *mcell1 = mcell_map[mcell]; // repalced by this one

    // save time
    std::set<int> times = cell_times_set_map[mcell];
    for (auto it1 = times.begin(); it1!=times.end(); it1++){
      cell_times_set_map[mcell1].insert(*it1);
      time_cells_set_map[*it1].insert(mcell1);
      time_cells_set_map[*it1].erase(mcell);
    }

    // Now remove ...
    cell_times_set_map.erase(mcell);
    auto it2 = find(mcells.begin(), mcells.end(), mcell);
    mcells.erase(it2);
    delete mcell;
  }
  to_be_removed.clear();
  
  //std::cout << mcells.size() << " " << to_be_saved.size() << " " << to_be_removed.size() << std::endl;
}


SMGCSelection PR3DCluster::Is_Connected(PR3DCluster* cluster1 , int offset){
  SMGCSelection temp_mcells;

  for (auto it = mcells.begin(); it!= mcells.end();it++){
    SlimMergeGeomCell *mcell = *it;
    int time_slice = mcell->GetTimeSlice();
    std::map<int,SMGCSet>& time_mcell_map = cluster1->get_time_cells_set_map();
    // std::cout << time_slice << " " << time_mcell_map.size() << std::endl;
    // std::vector<int> times;
    // if (time_cells_set_map.find(time_slice)!=time_cells_set_map.end())
    //   times.push_back(time_slice);
    // if (time_cells_set_map.find(time_slice-1)!=time_cells_set_map.end())
    //   times.push_back(time_slice-1);
    // if (time_cells_set_map.find(time_slice+1)!=time_cells_set_map.end())
    //   times.push_back(time_slice+1);

    bool flag = false;
    // for (int i=0;i!=times.size();i++){
    if (time_mcell_map.find(time_slice)!=time_mcell_map.end()){
      SMGCSet& mcells1 = time_mcell_map[time_slice];
      for (auto it1 = mcells1.begin(); it1!=mcells1.end();it1++){
   	SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)(*it1);
   	if (mcell->Overlap_fast(mcell1,offset)){
   	  flag = true;
   	  break;
   	}
      }
      //     //  if (flag) break;
    }
    if (flag) temp_mcells.push_back(mcell);
  }
  
  return temp_mcells;
}
