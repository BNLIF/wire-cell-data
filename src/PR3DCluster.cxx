#include "WireCellData/PR3DCluster.h"
#include "WireCellData/TPCParams.h"
#include "WireCellData/Singleton.h"

#include "TMatrixDEigen.h"
#include "TH2F.h"
#include "TVector3.h"

#include <boost/graph/connected_components.hpp>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>


using namespace WireCell;

PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
  point_cloud = 0;
  graph = 0;
  source_wcp_index = -1;
  flag_fine_tracking = false;
}

PR3DCluster::~PR3DCluster(){
  if (point_cloud!=(ToyPointCloud*)0)
    delete point_cloud;
  if (graph!=(MCUGraph*)0)
    delete graph;
}

std::vector<int> PR3DCluster::get_uvwt_range(){
  std::set<int> set_u, set_v, set_w, set_t;
  for (auto it=mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = *it;
    int time_slice = mcell->GetTimeSlice();
    set_t.insert(time_slice);
    GeomWireSelection uwires = mcell->get_uwires();
    GeomWireSelection vwires = mcell->get_vwires();
    GeomWireSelection wwires = mcell->get_wwires();
    for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
      const GeomWire* wire = (*it1);
      set_u.insert(wire->index());
    }
    for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
      const GeomWire* wire = (*it1);
      set_v.insert(wire->index());
    }
    for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
      const GeomWire* wire = (*it1);
      set_w.insert(wire->index());
    }
  }
  std::vector<int> results;
  results.push_back(set_u.size());
  results.push_back(set_v.size());
  results.push_back(set_w.size());
  results.push_back(set_t.size());
  
  return results;
}

std::pair<int,int> PR3DCluster::get_num_points(Point& p, TVector3& dir){
  int num_p1 = 0;
  int num_p2 = 0;

  // loop through all the points
  const int N = point_cloud->get_num_points();
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  for (int i=0;i!=N;i++){
    TVector3 dir1(cloud.pts[i].x - p.x, cloud.pts[i].y - p.y, cloud.pts[i].z - p.z);
    if (dir1.Dot(dir)>=0){
      num_p1++;
    }else{
      num_p2++;
    }
  }
  
  return std::make_pair(num_p1,num_p2);
}


void PR3DCluster::Update_mcell_cluster_map(std::map<WireCell::SlimMergeGeomCell*,WireCell::PR3DCluster*>& mcell_cluster_map){
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    mcell_cluster_map[mcell] = this;
  }
}
 
void PR3DCluster::Create_point_cloud(WireCell::ToyPointCloud *global_point_cloud){
  if (point_cloud!=(ToyPointCloud*)0)
    return;
  
  point_cloud = new ToyPointCloud();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    PointVector pts = mcell->get_sampling_points();

    if (global_point_cloud!=(ToyPointCloud*)0)
      global_point_cloud->AddPoints(pts,mcell->get_sampling_points_wires(),mcell);  
    
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
  WireCell::WC2DPointCloud<double>& cloud_u = point_cloud->get_cloud_u();
  WireCell::WC2DPointCloud<double>& cloud_v = point_cloud->get_cloud_v();
  WireCell::WC2DPointCloud<double>& cloud_w = point_cloud->get_cloud_w();
  

  
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
	pt_clouds.at(component[i])->AddPoint(cloud.pts[i],cloud_u.pts[i],cloud_v.pts[i],cloud_w.pts[i]);
	//   std::cout << "Vertex " << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " << cloud.pts[i].mcell << " " << cloud.pts[i].mcell->GetTimeSlice()  << " is in component " << component[i] << std::endl;
      }
      for (int j=0;j!=num;j++){
	pt_clouds.at(j)->build_kdtree_index();
      }
      // connect these graphs according to closest distance some how ...
      std::tuple<int,int,double> index_index_dis[num][num];
      for (int j=0;j!=num;j++){
	for (int k=j+1;k!=num;k++){
	  index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(pt_clouds.at(k));
	  //std::cout << j << " "<< k << std::get<0>(index_index_dis[j][k]) << " " << std::get<1>(index_index_dis[j][k]) << " " << std::get<2>(index_index_dis[j][k]) << std::endl;
	}
      }

      for (int j=0;j!=num-1;j++){
	// find the minimum
	std::tuple<int,int,double> min_dis = index_index_dis[j][j+1];
	int min_index = j+1;
	for (int k=j+2;k<num;k++){
	  if (std::get<2>(index_index_dis[j][k]) < std::get<2>(min_dis)){
	    min_dis = index_index_dis[j][k];
	    min_index = k;
	  }
	}

	for (int k=j+1;k<min_index;k++){
	  if (std::get<2>(index_index_dis[j][k]) < std::get<2>(index_index_dis[k][min_index]) ){
	    index_index_dis[k][min_index] = index_index_dis[j][k];
	  }
	}
	for (int k=min_index+1;k<num;k++){
	  if (std::get<2>(index_index_dis[j][k]) < std::get<2>(index_index_dis[min_index][k]) ){
	    index_index_dis[min_index][k]=index_index_dis[j][k];
	  }
	}
	//std::cout << j << " " << min_index << " " << std::get<0>(min_dis) << " " << std::get<1>(min_dis) << " "<< std::get<2>(min_dis) << std::endl;

	auto edge = add_edge(std::get<0>(min_dis),std::get<1>(min_dis),*graph);
	if (edge.second){
	  (*graph)[edge.first].dist = std::get<2>(min_dis);
	}
      }
      
      for (int i=0;i!=num;i++){
	delete pt_clouds.at(i);
      }
    }
  }

  // {
  //   std::vector<int> component(num_vertices(*graph));
  //   const int num = connected_components(*graph,&component[0]);
  //   if (num>1) std::cout << "Wrong! " << num << std::endl;
  // }
  
}

void PR3DCluster::dijkstra_shortest_paths(WCPointCloud<double>::WCPoint& wcp){
  if (graph==(MCUGraph*)0)
    Create_graph();
  if (wcp.index==source_wcp_index)
    return;
  source_wcp_index = wcp.index;
  parents.resize(num_vertices(*graph));
  distances.resize(num_vertices(*graph));
  
  auto v0 = vertex(wcp.index,*graph);
  boost::dijkstra_shortest_paths(*graph, v0,
				 weight_map(get(&EdgeProp::dist, *graph))
				 .predecessor_map(&parents[0])
				 .distance_map(&distances[0])
				 );
}


void PR3DCluster::cal_shortest_path(WCPointCloud<double>::WCPoint& wcp_target){
  dest_wcp_index = wcp_target.index;
  path_wcps.clear();
  path_mcells.clear();

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  int prev_i = -1;
  for(int i = dest_wcp_index; i!=source_wcp_index; i = parents[i]) {
    if (path_wcps.size()==0){
      path_wcps.push_front(cloud.pts[i]);
      path_mcells.push_front(cloud.pts[i].mcell);
    }else{
      path_wcps.push_front(cloud.pts[i]);
      if (cloud.pts[i].mcell!=path_mcells.front())
	path_mcells.push_front(cloud.pts[i].mcell);
    }
    if (i==prev_i) break;
    prev_i = i;
  }
  path_wcps.push_front(cloud.pts[source_wcp_index]);
  if (cloud.pts[source_wcp_index].mcell!=path_mcells.front())
    path_mcells.push_front(cloud.pts[source_wcp_index].mcell);
  
}

void PR3DCluster::fine_tracking(int num_pts_cut){
  // cut ... 
  if (path_wcps.size() < num_pts_cut) return;

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  
  //skip anypoint which is further away than 0.5 cm
  double low_dis_limit = 0.5*units::cm;
  std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec;
  std::vector<double> distances;
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      path_wcps_vec.push_back(*it);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	path_wcps_vec.push_back(*it);
	distances.push_back(dis);
      }
    }
    //    if (path_wcps_vec.size()==2) break;
  }
  // for (size_t i=0;i!=distances.size();i++){
  //   std::cout << i << " " << distances.at(i)/2./units::cm << std::endl;
  // }

  //form a map, (U,T) --> charge and error
  //form a map, (V,T) --> charge and error
  //form a map, (Z,T) --> charge and error
  std::map<std::pair<int,int>,double> map_2D_ut_charge, map_2D_ut_charge_err;
  std::map<std::pair<int,int>,double> map_2D_vt_charge, map_2D_vt_charge_err;
  std::map<std::pair<int,int>,double> map_2D_wt_charge, map_2D_wt_charge_err;
  for (auto it=mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    WireChargeMap& wire_charge_err_map = mcell->get_wirechargeerr_map();
    for (auto it = wire_charge_map.begin(); it!= wire_charge_map.end(); it++){
      const GeomWire* wire = it->first;
      double charge = it->second;
      double charge_err = wire_charge_err_map[wire];
      //std::cout << wire << " " << charge << " " << charge_err << std::endl;

      // hack the charge ... 
      if (charge <=0){
	continue;
	//charge = 1000;
	//charge_err = 1000;
      }
      
      if (wire->iplane()==0){
	map_2D_ut_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_ut_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }else if (wire->iplane()==1){
	map_2D_vt_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_vt_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }else{
	map_2D_wt_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_wt_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }
	
    }
    //    std::cout << wire_charge_map.size() << " " << wire_charge_err_map.size() << std::endl;
  }


  // Now prepare all the maps ...
  // map 3D index to set of 2D points
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DU_set;
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DV_set;
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DW_set;
  // map 2D points to 3D indices
  std::map<std::pair<int,int>,std::set<int>> map_2DU_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DV_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DW_3D_set;
  
  //std::cout << map_2DU_index.size() << " " << map_2DV_index.size() << " " << map_2DW_index.size() << std::endl;

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  
  float coef1 = 2 * pow(sin(angle_u),2);
  float coef2 = 2 * (pow(sin(angle_u),2) - pow(cos(angle_u),2));
  
  //std::cout << pitch_u/units::cm << " " << pitch_v/units::cm << " " << pitch_w/units::cm << " " << time_slice_width/units::cm << std::endl;
  // Loop any point and try to find its 3-level neibours ...
  typedef boost::property_map<MCUGraph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index,*graph);
  typedef boost::graph_traits<MCUGraph>::adjacency_iterator adjacency_iterator;
  int nlevel = 3;
  double dis_cut,time_cut;
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    int current_index = path_wcps_vec.at(i).index;
    if (i==0){
      dis_cut = std::min(distances.at(i) * 0.6,0.8*units::cm);
    }else if (i==path_wcps_vec.size()-1){
      dis_cut = std::min(distances.back() * 0.6,0.8*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*0.9,distances.at(i)*0.9),1.2*units::cm);
    }
    
    time_cut = 3; // allow +- 3 time slices and then distance cut ... 
    std::set<std::pair<int,int>> T2DU_set;
    std::set<std::pair<int,int>> T2DV_set;
    std::set<std::pair<int,int>> T2DW_set;
    map_3D_2DU_set[i] = T2DU_set;
    map_3D_2DV_set[i] = T2DV_set;
    map_3D_2DW_set[i] = T2DW_set;
    
    
    
    std::set<int> total_vertices_found;
    std::set<int> vertices_to_be_examined;
    std::set<int> vertices_saved_for_next;
    total_vertices_found.insert(current_index);
    vertices_to_be_examined.insert(current_index);

    for (int j=0;j!=nlevel;j++){
      for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
     	int temp_current_index = (*it);
     	std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph),*graph);
     	for (; neighbors.first!=neighbors.second; ++neighbors.first){
    	  //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
    	  if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
	    total_vertices_found.insert(index(*neighbors.first));
	    vertices_saved_for_next.insert(index(*neighbors.first));
    	  }
    	}
      }
      vertices_to_be_examined = vertices_saved_for_next;
    }
    SMGCSet nearby_mcells_set;
    for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
      SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
      nearby_mcells_set.insert(mcell);
    }
    // std::cout << i << " " << total_vertices_found.size() << " " << nearby_mcells_set.size() << std::endl;
    

    int cur_time_slice = cloud.pts[current_index].mcell->GetTimeSlice();
    int cur_wire_u = cloud.pts[current_index].index_u;
    int cur_wire_v = cloud.pts[current_index].index_v;
    int cur_wire_w = cloud.pts[current_index].index_w;

    // if (abs(cur_time_slice-1261)==0)
    //   std::cout << "Center: " << i << " " << path_wcps_vec.at(i).mcell << " " << current_index <<  " " << cloud.pts[current_index].index << " " <<
    // 	cloud.pts[current_index].mcell << " " << cloud.pts[current_index].x << " " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_v << std::endl;

    // Now fill the other maps ...
    for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int this_time_slice = mcell->GetTimeSlice();
      double rem_dis_cut = pow(dis_cut,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      if (rem_dis_cut >0 && fabs(cur_time_slice-this_time_slice)<=time_cut){
    	//	std::cout << this_time_slice << std::endl;
    	GeomWireSelection uwires = mcell->get_uwires();
    	GeomWireSelection vwires = mcell->get_vwires();
    	GeomWireSelection wwires = mcell->get_wwires();
	
    	float min_u_dis;
    	if (cur_wire_u < uwires.front()->index()){
    	  min_u_dis = uwires.front()->index()-cur_wire_u;
    	}else if (cur_wire_u >= uwires.front()->index() &&
    		  cur_wire_u <= uwires.back()->index()){
    	  min_u_dis = 0;
    	}else{
    	  min_u_dis = cur_wire_u-uwires.back()->index();
    	}
    	float min_v_dis;
    	if (cur_wire_v < vwires.front()->index()){
    	  min_v_dis = vwires.front()->index()-cur_wire_v;
    	}else if (cur_wire_v >= vwires.front()->index() &&
    		  cur_wire_v <= vwires.back()->index()){
    	  min_v_dis = 0;
    	}else{
    	  min_v_dis = cur_wire_v-vwires.back()->index();
    	}
    	float min_w_dis;
    	if (cur_wire_w < wwires.front()->index()){
    	  min_w_dis = wwires.front()->index()-cur_wire_w;
    	}else if (cur_wire_w >= wwires.front()->index() &&
    		  cur_wire_w <= wwires.back()->index()){
    	  min_w_dis = 0;
    	}else{
    	  min_w_dis = cur_wire_w-wwires.back()->index();
    	}

	// Note this assumes w is vertical wires, U and V have equal angle between then ...

    	float range_u = rem_dis_cut*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_v = rem_dis_cut*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_w = (rem_dis_cut*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	//std::cout << coef1 << " " << coef2 << std::endl;

    	// if (abs(cur_time_slice-1261)==0)
    	//   std::cout << min_u_dis << " " << min_v_dis << " " << min_w_dis << std::endl;
	
    	if ( range_u > 0 && range_v >0 && range_w > 0){
    	  float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
    	  float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
    	  float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
    	  float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
    	  float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
    	  float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;

    	  // if (abs(cur_time_slice-1261)==0)
    	  //   std::cout << low_u_limit << " " << high_u_limit << " " << low_v_limit << " " << high_v_limit << " " << low_w_limit << " " << high_w_limit << std::endl;
	  
    	  WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    	  for (auto it1 = wire_charge_map.begin(); it1!= wire_charge_map.end(); it1++){
	    const GeomWire *wire = it1->first;
	    if (it1->second >0){
	      //if (1>0){
	      if (wire->iplane()==0){
		// U plane ...
		if (wire->index() >= low_u_limit && wire->index() <= high_u_limit){
		  map_3D_2DU_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DU_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DU_3D_set.end()){
		    std::set<int>  temp_set;
		    temp_set.insert(i);
		    map_2DU_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
		    map_2DU_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }else if (wire->iplane()==1){
		// V plane ...
		if (wire->index() >= low_v_limit && wire->index() <= high_v_limit){
		  map_3D_2DV_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DV_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DV_3D_set.end()){
		    std::set<int>  temp_set;
		    temp_set.insert(i);
		    map_2DV_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
    		  map_2DV_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }else{
		// W plane ... 
		if (wire->index() >= low_w_limit && wire->index() <= high_w_limit){
		  map_3D_2DW_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DW_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DW_3D_set.end()){
		    std::set<int> temp_set;
		    temp_set.insert(i);
		    map_2DW_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
		    map_2DW_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }
	    }
    	  }
	}
      }
    }
  } // i loop ... 


  // map 2D points to its index
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DU_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DV_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DW_index;

  for (auto it = map_2DU_3D_set.begin(); it!= map_2DU_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DU_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DV_3D_set.begin(); it!= map_2DV_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DV_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DW_3D_set.begin(); it!= map_2DW_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DW_index.push_back(std::make_pair(it->first,*it1));
    }
  }

  //  std::cout << vec_2DU_index.size() << " " << vec_2DV_index.size() << " " << vec_2DW_index.size() << std::endl;
  

  // for (auto it = map_2D_ut_charge.begin(); it!=map_2D_ut_charge.end();it++){
  //   std::cout << "U_2: " << it->first.first << " " << it->first.second << std::endl;
  // }
  // for (auto it = map_2D_vt_charge.begin(); it!=map_2D_vt_charge.end();it++){
  //   std::cout << "V_2: " << it->first.first << " " << it->first.second << std::endl;
  // }
  // for (auto it = map_2D_wt_charge.begin(); it!=map_2D_wt_charge.end();it++){
  //   std::cout << "W_2: " << it->first.first << " " << it->first.second << std::endl;
  // }
  // for (size_t i=0;i!=path_wcps_vec.size();i++){
  //   std::cout << i << " " << path_wcps_vec.at(i).mcell->GetTimeSlice() << std::endl;
  //   std::cout << "U_3: " << path_wcps_vec.at(i).index_u << " " << path_wcps_vec.at(i).x << std::endl;
  //   std::cout << "V_3: " << path_wcps_vec.at(i).index_v << " " << path_wcps_vec.at(i).x << std::endl;
  //   std::cout << "W_3: " << path_wcps_vec.at(i).index_w << " " << path_wcps_vec.at(i).x << std::endl;
  // }
  // std::cout << map_2DU_index.size() << " " << map_2D_ut_charge.size() << " "
  // 	    << map_2DV_index.size() << " " << map_2D_vt_charge.size() << " "
  //   	    << map_2DW_index.size() << " " << map_2D_wt_charge.size() << " "
  // 	    << path_wcps_vec.size() << std::endl;
  

  
  // Now start to prepare the matrix and solve it ...
  // need to establish how to convert X position to T index
  // T = slope_x * ( x + offset_x);
  double slope_x = 1./time_slice_width;
  double first_t_dis = path_wcps_vec.at(0).mcell->GetTimeSlice()*time_slice_width - path_wcps_vec.at(0).x;
  double offset_t = first_t_dis / time_slice_width;
  
  
  
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w;
  double offset_u = -first_u_dis/pitch_u;
  double offset_v = -first_v_dis/pitch_v;
 
  //  for (size_t i=0;i!=path_wcps_vec.size();i++){
    //std::cout << i << " " << path_wcps_vec.at(i).mcell->GetTimeSlice() - offset_t << " " << slope_x * path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).mcell->GetTimeSlice() - offset_t - slope_x * path_wcps_vec.at(i).x << std::endl;
    //std::cout << i << " " << path_wcps_vec.at(i).index_w - offset_w << " " << slope_zw * path_wcps_vec.at(i).z << " " << path_wcps_vec.at(i).index_w - offset_w-slope_zw * path_wcps_vec.at(i).z <<  std::endl;
    //std::cout << i << " " << path_wcps_vec.at(i).index_u - offset_u << " " << slope_yu * path_wcps_vec.at(i).y + slope_zu * path_wcps_vec.at(i).z << " " << path_wcps_vec.at(i).index_u - offset_u - (slope_yu * path_wcps_vec.at(i).y + slope_zu * path_wcps_vec.at(i).z) << std::endl;
    // std::cout << i << " " << path_wcps_vec.at(i).index_v - offset_v << " " << slope_yv * path_wcps_vec.at(i).y + slope_zv * path_wcps_vec.at(i).z << " " << path_wcps_vec.at(i).index_v - offset_v - (slope_yv * path_wcps_vec.at(i).y + slope_zv * path_wcps_vec.at(i).z) << std::endl;
  //}

  
  // form matrix ...
  int n_3D_pos = 3 * path_wcps_vec.size();
  int n_2D_u = 2 * vec_2DU_index.size();
  int n_2D_v = 2 * vec_2DV_index.size();
  int n_2D_w = 2 * vec_2DW_index.size();
  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    pos_3D_init(3*i) = path_wcps_vec.at(i).x;
    pos_3D_init(3*i+1) = path_wcps_vec.at(i).y;
    pos_3D_init(3*i+2) = path_wcps_vec.at(i).z;
  }
  
  // fill in the measurement ...
  for (size_t index = 0; index!=vec_2DU_index.size(); index++){
    double charge = map_2D_ut_charge[vec_2DU_index.at(index).first];
    double charge_err = map_2D_ut_charge_err[vec_2DU_index.at(index).first];
    int n_divide = map_2DU_3D_set[vec_2DU_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_u_2D(2*index) =  scaling * (vec_2DU_index.at(index).first.first - offset_u);
    data_u_2D(2*index+1) = scaling * (vec_2DU_index.at(index).first.second - offset_t);

    int index_3D = vec_2DU_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RU.insert(2*index,3*index_3D+1) = scaling * slope_yu; // Y--> U
    RU.insert(2*index,3*index_3D+2) = scaling * slope_zu; // Z--> U
    RU.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
    //std::cout << index << " " << index_3D << std::endl;
  }
  for (size_t index = 0; index!=vec_2DV_index.size(); index++){
    double charge = map_2D_vt_charge[vec_2DV_index.at(index).first];
    double charge_err = map_2D_vt_charge_err[vec_2DV_index.at(index).first];
    int n_divide = map_2DV_3D_set[vec_2DV_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_v_2D(2*index) =  scaling * (vec_2DV_index.at(index).first.first - offset_v);
    data_v_2D(2*index+1) = scaling * (vec_2DV_index.at(index).first.second - offset_t);

    int index_3D = vec_2DV_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RV.insert(2*index,3*index_3D+1) = scaling * slope_yv; // Y--> V
    RV.insert(2*index,3*index_3D+2) = scaling * slope_zv; // Z--> V
    RV.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
    //std::cout << index << " " << index_3D << std::endl;
  }
  for (size_t index = 0; index!=vec_2DW_index.size(); index++){
    double charge = map_2D_wt_charge[vec_2DW_index.at(index).first];
    double charge_err = map_2D_wt_charge_err[vec_2DW_index.at(index).first];
    int n_divide = map_2DW_3D_set[vec_2DW_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_w_2D(2*index) =  scaling * (vec_2DW_index.at(index).first.first - offset_w);
    data_w_2D(2*index+1) = scaling * (vec_2DW_index.at(index).first.second - offset_t);

    int index_3D = vec_2DW_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RW.insert(2*index,3*index_3D+2) = scaling * slope_zw; // Z--> W
    RW.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
  }
 

 //  for (int k=0;k<RV.outerSize();++k){
 //    for (Eigen::SparseMatrix<double>::InnerIterator it(RV,k); it; ++it){
 //      std::cout << it.value() << " " << it.row() << " "<< it.col() << " " << it.index() << std::endl;
 //    }
 // }
  
  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());

  // double ave_distance = 0;
  // for (size_t i=0;i!=distances.size();i++){
  //   ave_distance += distances.at(i);
  // }
  // ave_distance /= distances.size();

  
  // get initial chi2 ...
  // Eigen::VectorXd chi_u = data_u_2D - RU * pos_3D_init;
  // Eigen::VectorXd chi_v = data_v_2D - RV * pos_3D_init;
  // Eigen::VectorXd chi_w = data_w_2D - RW * pos_3D_init;

  // double chi2 = chi_u.squaredNorm() + chi_v.squaredNorm() + chi_w.squaredNorm();
  // std::cout << sqrt(chi2 /(path_wcps_vec.size() * 1.)) << " " << lambda/units::cm*3 << std::endl;
  
  double lambda = 2 // strength 
    * sqrt(9. //  average chi2 guessted
 	   * (vec_2DU_index.size() + vec_2DV_index.size() + vec_2DW_index.size()) // how many of them
 	   * 6 * 6 // charge/charge_err estimation ... 
 	   /(path_wcps_vec.size() * 1.)); //weighting
  lambda = 0;
  //double dis_range = 0.5*units::cm/sqrt(3.);
  double angle_range = 0.25;
 
  //std::cout << lambda/ dis_range << std::endl;
  
  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos) ;
  //Eigen::SparseMatrix<double> PMatrix(n_3D_pos, n_3D_pos) ;
  // distances[i]
  // 2nd order ...
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    // PMatrix.insert(3*i,3*i)=1;
    // PMatrix.insert(3*i+1,3*i+1)=1;
    // PMatrix.insert(3*i+2,3*i+2)=1;
    if (i==0){
      FMatrix.insert(0,0) = -1./distances.at(0); // X
      FMatrix.insert(0,3) = 1./distances.at(0);
      FMatrix.insert(1,1) = -1./distances.at(0); // Y
      FMatrix.insert(1,4) = 1./distances.at(0);
      FMatrix.insert(2,2) = -1./distances.at(0); // Z
      FMatrix.insert(2,5) = 1./distances.at(0);
    }else if (i==path_wcps_vec.size()-1){
      FMatrix.insert(3*i,3*i) = -1./distances.at(path_wcps_vec.size()-2); // X
      FMatrix.insert(3*i,3*i-3) = 1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(path_wcps_vec.size()-2);
    }else{
      FMatrix.insert(3*i,3*i-3) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i,3*i) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i,3*i+3) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i+1,3*i+4) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i+2,3*i+5) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
    }
  }

   //0th order
  // for (size_t i=0;i!=path_wcps_vec.size();i++){
  //   if (i==0){
  //     FMatrix.insert(0,0) = -1;
  //     FMatrix.insert(0,3) = 1.;
  //     FMatrix.insert(1,1) = -1.;
  //     FMatrix.insert(1,4) = 1.;
  //     FMatrix.insert(2,2) = -1.;
  //     FMatrix.insert(2,5) = 1.;
  //   }else if (i==path_wcps_vec.size()-1){
  //     FMatrix.insert(3*i,3*i) = -1.;
  //     FMatrix.insert(3*i,3*i-3) = 1.;
  //     FMatrix.insert(3*i+1,3*i+1) = -1.;
  //     FMatrix.insert(3*i+1,3*i-2) = 1.;
  //     FMatrix.insert(3*i+2,3*i+2) = -1.;
  //     FMatrix.insert(3*i+2,3*i-1) = 1.;
  //   }else{
  //     FMatrix.insert(3*i,3*i) = -1.;
  //     FMatrix.insert(3*i,3*i+3) = 1.;
  //     FMatrix.insert(3*i+1,3*i+1) = -1.;
  //     FMatrix.insert(3*i+1,3*i+4) = 1.;
  //     FMatrix.insert(3*i+2,3*i+2) = -1.;
  //     FMatrix.insert(3*i+2,3*i+5) = 1.;
  //   }
  // }
  
  // 2nd order
  // for (size_t i=0;i!=path_wcps_vec.size();i++){
  //   if (i==0){
  //     FMatrix.insert(0,0) = -1;
  //     FMatrix.insert(0,3) = 1.;
  //     FMatrix.insert(1,1) = -1.;
  //     FMatrix.insert(1,4) = 1.;
  //     FMatrix.insert(2,2) = -1.;
  //     FMatrix.insert(2,5) = 1.;
  //   }else if (i==path_wcps_vec.size()-1){
  //     FMatrix.insert(3*i,3*i) = -1.;
  //     FMatrix.insert(3*i,3*i-3) = 1.;
  //     FMatrix.insert(3*i+1,3*i+1) = -1.;
  //     FMatrix.insert(3*i+1,3*i-2) = 1.;
  //     FMatrix.insert(3*i+2,3*i+2) = -1.;
  //     FMatrix.insert(3*i+2,3*i-1) = 1.;
  //   }else{
  //     FMatrix.insert(3*i,3*i-3) = 1.;
  //     FMatrix.insert(3*i,3*i) = -2.;
  //     FMatrix.insert(3*i,3*i+3) = 1.;
  //     FMatrix.insert(3*i+1,3*i-2) = 1.;
  //     FMatrix.insert(3*i+1,3*i+1) = -2.;
  //     FMatrix.insert(3*i+1,3*i+4) = 1.;
  //     FMatrix.insert(3*i+2,3*i-1) = 1.;
  //     FMatrix.insert(3*i+2,3*i+2) = -2.;
  //     FMatrix.insert(3*i+2,3*i+5) = 1.;
  //   }
  // }

  
  FMatrix *= lambda/angle_range ; // disable the angue cut ... 
  
  
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  // Eigen::SparseMatrix<double> PMatrixT = Eigen::SparseMatrix<double>(PMatrix.transpose());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;// + PMatrixT * pos_3D_init * pow(lambda/dis_range,2);
  Eigen::SparseMatrix<double> A =   RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;// + PMatrixT * PMatrix * pow(lambda/dis_range,2);
  solver.compute(A);
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);
  
  //std::cout << "#iterations: " << solver.iterations() << std::endl;
  //std::cout << "#estimated error: " << solver.error() << std::endl;
  //std::cout << path_wcps_vec.size() << " " << map_2DU_index.size() << " " << map_2DV_index.size() << " " << map_2DW_index.size() << std::endl;

  flag_fine_tracking = true;
  fine_tracking_path.clear();
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    Point p;
    p.x = pos_3D(3*i);
    p.y = pos_3D(3*i+1);
    p.z = pos_3D(3*i+2);
    fine_tracking_path.push_back(p);
    // std::cout << p.x << " " << p.y << " " << p.z << " " << path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).y << " " << path_wcps_vec.at(i).z << std::endl;
  }



  
  // Eigen::VectorXd data_v_2D_p1 = RVT * RV * pos_3D1;
  // Eigen::VectorXd data_v_2D_p2 = RVT * RV * pos_3D;
  // Eigen::VectorXd data_v_2D_p3 = RVT * data_v_2D;
  // for (size_t i=0;i!=path_wcps_vec.size();i++){
  //   std::cout << i << " " << data_v_2D_p3(3*i) << " " << data_v_2D_p1(3*i) << " " << data_v_2D_p2(3*i)
  // 	      << " " << data_v_2D_p3(3*i+1) << " " << data_v_2D_p1(3*i+1) << " " << data_v_2D_p2(3*i+1)
  // 	      << " " << data_v_2D_p3(3*i+2) << " " << data_v_2D_p1(3*i+2) << " " << data_v_2D_p2(3*i+2) << std::endl;
  // }
  // for (size_t i=0;i!= map_2DV_index.size();i++){
  //   std::cout << i << " X " << data_v_2D(2*i) << " " << data_v_2D_p1(2*i) << " " << data_v_2D_p2(2*i) << " "
  // 	      << data_v_2D(2*i+1) << " " << data_v_2D_p1(2*i+1) << " " << data_v_2D_p2(2*i+1) << " " << std::endl;
  // }
  // Eigen::VectorXd data_v_2D_p1 = RV * pos_3D1;
  // for (size_t i=0;i!=n_2D_v;i++){
  //   std::cout << data_v_2D(i) << " " << data_v_2D_p1(i) << std::endl;
  // }
  
  // for (int k=0;k<RV.outerSize();++k){
  //   for (Eigen::SparseMatrix<double>::InnerIterator it(RV,k); it; ++it){
  //     std::cout << it.value() << " " << it.row() << " "<< it.col() << " " << it.index() << std::endl;
  //   }
  // }

  // for (size_t i=0;i!=path_wcps_vec.size();i++){
  //   std::cout << pos_3D(3*i) << " " << pos_3D(3*i+1) << " " << pos_3D(3*i+2) << " " << path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).y << " " << path_wcps_vec.at(i).z << std::endl;
  // }
  
  
  // Eigen::VectorXd test1(3); test1(0)=1;test1(1)=2; test1(2)=3;
  // Eigen::MatrixXd test2(2,3);
  // test2(0,0)=1; 
  // test2(1,1) = 1; test2(1,2) = 1;
  // Eigen::VectorXd test3 = test2 * test1;
  // std::cout << test3(0) << " " << test3(1) << std::endl;
  
  
  // copy the list into a vector...
  //std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec(std::begin(path_wcps), std::end(path_wcps));
  
  
  
  
  // go along the vector, and create a new set of points, according to prechosen Dt
  //double dt = 0.5*units::cm;

  // max of( a certain distance (predefined) and three times the closest distance)
  // save the corresponding mcells for these points as candidate 
  
  // current point, nearby points, and points satisfying a certain distance cut as candidates
  



  
  
  
  
  // need a Kalman filter like alg to organize the points???
  // std::vector<bool> flag_prev_dirs, flag_next_dirs;
  // flag_prev_dirs.resize(path_wcps_vec.size(), true);
  // flag_next_dirs.resize(path_wcps_vec.size(), true);

  // TVector3 prev_pos,next_pos,dir,dir1;
  // for (size_t j=0; j!=path_wcps_vec.size(); j++){
  //   int current_index = j;
    
  //   // calculate the direction vector
  //   prev_pos.SetXYZ(0,0,0);
  //   int num_pos = 0;
  //   for (int i=std::max(0,current_index - 10);i!=current_index+1;i++){
  //     prev_pos.SetXYZ(prev_pos.X() - path_wcps_vec.at(i).x,
  // 		      prev_pos.Y() - path_wcps_vec.at(i).y,
  // 		      prev_pos.Z() - path_wcps_vec.at(i).z);
  //     num_pos ++;
  //   }
  //   prev_pos *= 1./num_pos;

  //   next_pos.SetXYZ(0,0,0);
  //   num_pos = 0;
  //   for (int i=current_index; i!=std::min(current_index+11,int(path_wcps_vec.size()));i++){
  //     next_pos.SetXYZ(next_pos.X() + path_wcps_vec.at(i).x,
  // 		      next_pos.Y() + path_wcps_vec.at(i).y,
  // 		      next_pos.Z() + path_wcps_vec.at(i).z
  // 		      );
  //     num_pos++;
  //   }
  //   next_pos *= 1./num_pos;
  //   dir = next_pos - prev_pos;

  //   if (j>0){
  //     dir1.SetXYZ(path_wcps_vec.at(j).x - path_wcps_vec.at(j-1).x,
  // 		  path_wcps_vec.at(j).y - path_wcps_vec.at(j-1).y,
  // 		  path_wcps_vec.at(j).z - path_wcps_vec.at(j-1).z);
  //     if (dir1.Dot(dir)<0) flag_prev_dirs.at(j) = false;
  //   }
  //   if (j+1<path_wcps_vec.size()){
  //     dir1.SetXYZ(path_wcps_vec.at(j+1).x - path_wcps_vec.at(j).x,
  // 		  path_wcps_vec.at(j+1).y - path_wcps_vec.at(j).y,
  // 		  path_wcps_vec.at(j+1).z - path_wcps_vec.at(j).z);
  //     if (dir1.Dot(dir)<0) flag_next_dirs.at(j) = false;
  //   }
    
  // }

  // for (size_t i=0; i!=flag_prev_dirs.size(); i++){
  //   std::cout << i << " " << flag_prev_dirs.at(i) << " " << flag_next_dirs.at(i) << " " <<
  //     path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).y << " " << path_wcps_vec.at(i).z << std::endl;
  // }
  
   
}

TVector3 PR3DCluster::get_ft_dir_end(float mean_dis, float dis_cut){
  TVector3 dir(0,0,0);
  for (size_t i=1; i<fine_tracking_path.size();i++){
    float dis = sqrt(pow(fine_tracking_path.at(i).x - fine_tracking_path.at(0).x,2) + pow(fine_tracking_path.at(i).y - fine_tracking_path.at(0).y,2) + pow(fine_tracking_path.at(i).z - fine_tracking_path.at(0).z,2));
    if ( dis < dis_cut ){
      dir.SetX(dir.X() + (fine_tracking_path.at(i).x - fine_tracking_path.at(0).x) * exp(-dis/mean_dis));
      dir.SetY(dir.Y() + (fine_tracking_path.at(i).y - fine_tracking_path.at(0).y) * exp(-dis/mean_dis));
      dir.SetZ(dir.Z() + (fine_tracking_path.at(i).z - fine_tracking_path.at(0).z) * exp(-dis/mean_dis));
    }else{
      break;
    }
  }
  
  return dir.Unit();
}


std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> PR3DCluster::get_highest_lowest_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].y > highest_wcp.y)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].y < lowest_wcp.y)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(highest_wcp,lowest_wcp);
}

std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> PR3DCluster::get_front_back_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].z > highest_wcp.z)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].z < lowest_wcp.z)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(highest_wcp,lowest_wcp);
}


std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> PR3DCluster::get_earliest_latest_wcps(){
  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  WCPointCloud<double>::WCPoint highest_wcp = cloud.pts[0];
  WCPointCloud<double>::WCPoint lowest_wcp = cloud.pts[0];
  for (size_t i=1;i<cloud.pts.size();i++){
    if (cloud.pts[i].x > highest_wcp.x)
      highest_wcp = cloud.pts[i];
    if (cloud.pts[i].x < lowest_wcp.x)
      lowest_wcp = cloud.pts[i];
  }
  return std::make_pair(lowest_wcp,highest_wcp);
}


Point PR3DCluster::calc_ave_pos(Point&p, int N){
  std::map<WireCell::SlimMergeGeomCell*, Point> pts = point_cloud->get_closest_mcell(p,N);
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

std::pair<Point, double> PR3DCluster::get_closest_point_along_vec(Point& p_test1, TVector3 dir, double test_dis, double dis_step, double angle_cut, double dis_cut){
  bool flag = false;
  Point p_test;
  double min_dis = 1e9;
  double min_dis1 = 1e9;
  Point min_point = p_test1;
  for (int i=0; i!= int(test_dis/dis_step)+1;i++){
    p_test.x = p_test1.x + dir.X() * i * dis_step;
    p_test.y = p_test1.y + dir.Y() * i * dis_step;
    p_test.z = p_test1.z + dir.Z() * i * dis_step;
    
    std::pair<SlimMergeGeomCell*, Point> pts = get_closest_point_mcell(p_test);
    double dis = sqrt(pow(p_test.x - pts.second.x,2)+pow(p_test.y - pts.second.y,2)+pow(p_test.z - pts.second.z,2));
    double dis1 = sqrt(pow(p_test1.x - pts.second.x,2)+pow(p_test1.y - pts.second.y,2)+pow(p_test1.z - pts.second.z,2));
    if (dis < std::min(dis1 * tan(angle_cut/180.*3.1415926),dis_cut)){
      if (dis < min_dis){
	min_dis = dis;
	min_point = pts.second;
	min_dis1 = dis1;
      }
      if (dis < 3*units::cm)
	return std::make_pair(pts.second,dis1);
    }
  }

  return std::make_pair(min_point,min_dis1);
  
}

std::pair<SlimMergeGeomCell*, Point> PR3DCluster::get_closest_point_mcell(Point& p_test){
  std::map<SlimMergeGeomCell*,Point> pts = point_cloud->get_closest_mcell(p_test,1);
  return std::make_pair((*pts.begin()).first,(*pts.begin()).second);
}

int PR3DCluster::get_num_points(Point& p_test, double dis){
  return point_cloud->get_closest_points(p_test, dis).size();
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

void PR3DCluster::Calc_PCA(PointVector& points){
  center.x=0; center.y=0; center.z=0;
  int nsum = 0;
  for (auto it = mcells.begin(); it!=mcells.end();it++){
    for (int k=0;k!=points.size();k++){
      center.x += points.at(k).x;
      center.y += points.at(k).y;
      center.z += points.at(k).z;
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
	PointVector& ps = points;
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

  PCA_values[0] = eigen_values(0,0) ;
  PCA_values[1] = eigen_values(1,1) ;
  PCA_values[2] = eigen_values(2,2) ;
  
  // std::cout << eigen_values(0,0) << " " << eigen_values(1,1) << " " << eigen_values(2,2) << std::endl;
  //std::cout << eigen_vectors(0,0) << " " << eigen_vectors(0,1) << " " << eigen_vectors(0,2) << std::endl;
  for (int i=0;i!=3;i++){
    PCA_axis[i].x = eigen_vectors(0,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].y = eigen_vectors(1,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].z = eigen_vectors(2,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));

  }
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

  PCA_values[0] = eigen_values(0,0) ;
  PCA_values[1] = eigen_values(1,1) ;
  PCA_values[2] = eigen_values(2,2) ;
  
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
