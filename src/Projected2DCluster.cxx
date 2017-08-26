#include "WireCellData/Projected2DCluster.h"

using namespace WireCell;

Projected2DCluster::Projected2DCluster(WirePlaneType_t plane_no)
  : plane_no(plane_no)
{
  time_slice_limit[0] = -1;
  time_slice_limit[1] = -1;

  wire_limit[0] = -1;
  wire_limit[1] = -1;
}

Projected2DCluster::~Projected2DCluster(){
}

void Projected2DCluster::Print(){
  for (auto it = time_slice_array.begin(); it!= time_slice_array.end(); it++){
    std::cout << *it << ": ";
    for (auto it1 = time_wires_map[*it].begin(); it1!=time_wires_map[*it].end(); it1++){
      std::cout << "(" << it1->first << ", " << it1->second << ") ";
    }
    std::cout << std::endl;
  }
}


void Projected2DCluster::AddCell(SlimMergeGeomCell *mcell){

  std::vector<WirePlaneType_t> bad_planes = mcell->get_bad_planes();
  if (find(bad_planes.begin(),bad_planes.end(),plane_no)==bad_planes.end()){
    int time_slice = mcell->GetTimeSlice();

    GeomWireSelection wires;
    if (plane_no == WirePlaneType_t(0)){
      wires = mcell->get_uwires();
    }else if (plane_no == WirePlaneType_t(1)){
      wires = mcell->get_vwires();
    }else if (plane_no == WirePlaneType_t(2)){
      wires = mcell->get_wwires();
    }
    int low_wire = wires.front()->index();
    int high_wire = wires.back()->index();

    //std::cout << low_wire << " " << high_wire << std::endl;
    
    // update the limits .... 
    if (time_slice_limit[0]==-1){
      time_slice_limit[0] = time_slice;
      time_slice_limit[1] = time_slice;
      
      wire_limit[0] = low_wire;
      wire_limit[1] = high_wire;
    }else{
      if (time_slice < time_slice_limit[0])
	time_slice_limit[0] = time_slice;
      if (time_slice > time_slice_limit[1])
	time_slice_limit[1] = time_slice;
      if (low_wire < wire_limit[0])
	wire_limit[0] = low_wire;
      if (high_wire > wire_limit[1])
	wire_limit[1] = high_wire;
    }


    // update the contents
    auto it = find(time_slice_array.begin(), time_slice_array.end(), time_slice);
    if (it==time_slice_array.end()){
      // fill in for the first time ... 
      time_slice_array.push_back(time_slice);
      std::list<std::pair<int,int>> wire_pairs;
      wire_pairs.push_back(std::make_pair(low_wire,high_wire));
      time_wires_map[time_slice] = wire_pairs;
    }else{
      
      // no need to add anything in time_slice
      std::list<std::pair<int,int>>& wire_pairs = time_wires_map[time_slice];

      // need to add wires ...
      if (high_wire < wire_pairs.front().first){
	wire_pairs.push_front(std::make_pair(low_wire,high_wire));
      }else if (low_wire > wire_pairs.back().second){
	wire_pairs.push_back(std::make_pair(low_wire,high_wire));
      }else{
	// loop all of them

	std::list<std::pair<int,int>>::iterator insert_pos = wire_pairs.begin();
	std::vector<std::list<std::pair<int,int>>::iterator> to_be_removed;
	
	for (auto it=wire_pairs.begin(); it!=wire_pairs.end(); it++){
	  //judge overlapping
	  if (high_wire >= it->first && low_wire <= it->second){
	    to_be_removed.push_back(it);
	  }
	  
	  if (it->second < low_wire){
	    insert_pos = it;
	    insert_pos ++;
	  }
	  // auto it1 = it; // current element
	  // it1++; // next element
	}
	
	// update the low and high limit
	int temp_low_wire = low_wire;
	int temp_high_wire = high_wire;
	if (to_be_removed.size()>0){
	  if (to_be_removed.front()->first < temp_low_wire){
	    temp_low_wire = to_be_removed.front()->first;
	  }
	  if (to_be_removed.back()->second > temp_high_wire){
	    temp_high_wire = to_be_removed.back()->second;
	  }
	}
	wire_pairs.insert(insert_pos,std::make_pair(temp_low_wire,temp_high_wire));
	if (to_be_removed.size()==1){
	  wire_pairs.erase(to_be_removed.front());
	}else if (to_be_removed.size()>1){
	  auto it1 =  to_be_removed.back();
	  it1++;
	  wire_pairs.erase(to_be_removed.front(),it1);
	}
      }
    }
  }
}



int Projected2DCluster::judge_coverage(Projected2DCluster *cluster){
  // +1 cluster belong to this
  // -1 this belong to cluster
  // 0 they do not belong to each other
  // +2 they are identical with contents ...
  // -2 they are identical and both empty ...
  
  if (plane_no!=cluster->GetPlaneNo()){
    std::cout <<"Something Wrong! " << std::endl;
    return 0;
  }
  
  // if the cluster is empty ...
  if (cluster->get_number_time_slices()==0 && get_number_time_slices()==0){
    return -2; // all bad ... 
  }else if (cluster->get_number_time_slices()==0 && get_number_time_slices()!=0){
    return 1; // cluster is part of this
  }else if (cluster->get_number_time_slices()!=0 && get_number_time_slices()==0){
    return -1; // this is part of cluster
  }else{ // both are not empty
    bool is_this_inside_cluster = true;
    bool is_cluster_inside_this = true; 

    std::vector<int>& cluster_time_slice_array = cluster->get_time_slice_array();
    std::map<int,std::list<std::pair<int,int>>>& cluster_time_wires_map = cluster->get_time_wires_map();
    
    // loop content of this
    for (auto it = time_slice_array.begin(); it!= time_slice_array.end(); it++){
      int curr_time_slice = *it;
      std::list<std::pair<int,int>> wire_pairs = time_wires_map[curr_time_slice];
      if (cluster_time_wires_map.find(curr_time_slice)==cluster_time_wires_map.end()){
	is_this_inside_cluster = false;
	break;
      }else{
	// go into wire pairs ...
	for (auto it1 = wire_pairs.begin(); it1!= wire_pairs.end(); it1++){
	  int low_wire_limit = it1->first;
	  int high_wire_limit = it1->second;
	  is_this_inside_cluster=false;

	  for (auto it2 = cluster_time_wires_map[curr_time_slice].begin(); it2!= cluster_time_wires_map[curr_time_slice].end(); it2++){
	    int cluster_low_wire_limit = it2->first;
	    int cluster_high_wire_limit = it2->second;

	    if (low_wire_limit >= cluster_low_wire_limit &&
		high_wire_limit <= cluster_high_wire_limit){
	      is_this_inside_cluster = true;
	      break;
	    }else if (high_wire_limit>= cluster_low_wire_limit &&
		      low_wire_limit <= cluster_high_wire_limit){
	      break;
	    }
	  }
	  if (!is_this_inside_cluster)
	    break;
	}
      }
      if (!is_this_inside_cluster)
	break;
    }
    
    // loop content of the cluster
    for (auto it = cluster_time_slice_array.begin(); it!= cluster_time_slice_array.end(); it++){
      int curr_time_slice = *it;
      std::list<std::pair<int,int>> wire_pairs = cluster_time_wires_map[curr_time_slice];
      if (time_wires_map.find(curr_time_slice)==time_wires_map.end()){
	is_cluster_inside_this = false;
	break;
      }else{
	for (auto it1 = wire_pairs.begin(); it1!= wire_pairs.end(); it1++){
	  int cluster_low_wire_limit = it1->first;
	  int cluster_high_wire_limit = it1->second;
	  is_cluster_inside_this=false;


	  for (auto it2 = time_wires_map[curr_time_slice].begin(); it2!= time_wires_map[curr_time_slice].end(); it2++){
	    int low_wire_limit = it2->first;
	    int high_wire_limit = it2->second;

	    if (low_wire_limit <= cluster_low_wire_limit &&
	   	high_wire_limit >= cluster_high_wire_limit){
	      is_cluster_inside_this = true;
	      break;
	    }else if (cluster_high_wire_limit>= low_wire_limit &&
	   	      cluster_low_wire_limit <= high_wire_limit){
	      break;
	    }
	  }
	  if (!is_cluster_inside_this)
	    break;
	}
      }
      if (!is_cluster_inside_this)
	break;
    }
    
    
    if (is_this_inside_cluster && is_cluster_inside_this){
      return 2;
    }else if (is_this_inside_cluster && !is_cluster_inside_this){
      return -1;
    }else if (!is_this_inside_cluster && is_cluster_inside_this){
      return 1;
    }else{
      return 0;
    }
  }
  
  return 0;
}
