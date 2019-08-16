#include "WireCellData/Slim3DCluster.h"

using namespace std;
using namespace WireCell;

Slim3DCluster::Slim3DCluster(SlimMergeGeomCell &cell)
  :u_proj(0)
  ,v_proj(0)
  ,w_proj(0)
  ,flag_saved(0)
  ,flag_saved_1(0)
  ,flag_saved_u(0)
  ,flag_saved_v(0)
  ,flag_saved_w(0)
  ,id(-1)
  , total_charge (0)
  , min_total_charge(0)
{
  std::set<SlimMergeGeomCell*>  abc;
  abc.insert(&cell);
  cluster.push_back(abc);
}

Slim3DCluster::~Slim3DCluster(){
  cluster.clear();
  gcluster.clear();

  // front_cell_map.clear();
  // back_cell_map.clear();
  
  if (u_proj!=0) delete u_proj;
  if (v_proj!=0) delete v_proj;
  if (w_proj!=0) delete w_proj;
}

GeomCellSelection& Slim3DCluster::get_allcell(){
  gcluster.clear();
  for (int i=0;i!=cluster.size();i++){
    std::set<SlimMergeGeomCell*>  abc = cluster.at(i);
    for (auto it = abc.begin(); it!=abc.end(); ++ it){
      gcluster.push_back(*it);
    }
  }
  return gcluster;
}


GeomCellSelection Slim3DCluster::Is_Connected(Slim3DDeadCluster* cluster1 , int offset){
  GeomCellSelection mcells;
  for (auto it = gcluster.begin(); it!= gcluster.end();it++){
    SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
    int time_slice = mcell->GetTimeSlice();

    std::map<int,GeomCellSetp>& time_mcell_map = cluster1->get_cluster();
    // std::vector<int> times;
    // if (time_mcell_map.find(time_slice)!=time_mcell_map.end())
    //   times.push_back(time_slice);
    // if (time_mcell_map.find(time_slice-1)!=time_mcell_map.end())
    //   times.push_back(time_slice-1);
    // if (time_mcell_map.find(time_slice+1)!=time_mcell_map.end())
    //   times.push_back(time_slice+1);

    bool flag = false;
    
    // for (int i=0;i!=times.size();i++){
    if (time_mcell_map.find(time_slice)!=time_mcell_map.end()){
      GeomCellSetp& mcells1 = time_mcell_map[time_slice];
      for (auto it1 = mcells1.begin(); it1!=mcells1.end();it1++){
	SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)(*it1);
	if (mcell->Overlap_fast(mcell1,offset)){
	  flag = true;
	  break;
	}
      }
      //  if (flag) break;
    }
    if (flag) mcells.push_back(mcell);
  }
  
  return mcells;
}


Projected2DCluster* Slim3DCluster::get_projection(WirePlaneType_t plane){
  if (plane==WirePlaneType_t(0)){
    return u_proj;
  }else if(plane==WirePlaneType_t(1)){
    return v_proj;
  }else if(plane==WirePlaneType_t(2)){
    return w_proj;
  }
}

void Slim3DCluster::Calc_Projection(){
  u_proj = new Projected2DCluster(WirePlaneType_t(0),id);
  v_proj = new Projected2DCluster(WirePlaneType_t(1),id);
  w_proj = new Projected2DCluster(WirePlaneType_t(2),id);

  total_charge = 0;
  min_total_charge = 0;
  for (auto it = cluster.begin(); it!= cluster.end(); it++){
    for (auto it1 = (*it).begin(); it1!= (*it).end(); it1++){
      SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it1);
      total_charge += mcell->Estimate_total_charge();
      min_total_charge += mcell->Estimate_minimum_charge();
      u_proj->AddCell(mcell);
      v_proj->AddCell(mcell);
      w_proj->AddCell(mcell);
    }
  }

  u_proj->set_estimated_total_charge(min_total_charge);
  v_proj->set_estimated_total_charge(min_total_charge);
  w_proj->set_estimated_total_charge(min_total_charge);

  // std::cout << "Xin: U " << std::endl;
  // u_proj->Print();

  // std::cout << "Xin: V " << std::endl;
  // v_proj->Print();

  // std::cout << "Xin: W " << std::endl;
  // w_proj->Print();
}

void Slim3DCluster::DirectAddCell(SlimMergeGeomCell &cell){
  gcluster.push_back(&cell);
}
void Slim3DCluster::DirectOrderCells(){
  std::map<int,std::set<SlimMergeGeomCell*> > temp_map;
  for (auto it = gcluster.begin(); it!=gcluster.end(); it++){
    SlimMergeGeomCell* mcell = (SlimMergeGeomCell*)(*it);
    int time_slice = mcell->GetTimeSlice();
    if (temp_map.find(time_slice)==temp_map.end()){
      std::set<SlimMergeGeomCell*> temp_set;
      temp_set.insert(mcell);
      temp_map[time_slice] = temp_set;
    }else{
      temp_map[time_slice].insert(mcell);
    }
  }

  for (auto it = temp_map.begin(); it!=temp_map.end(); it++){
    cluster.push_back(it->second);
  }
  
}

int Slim3DCluster::AddCell(SlimMergeGeomCell &cell,int offset){
  std::set<SlimMergeGeomCell*>&  abc = cluster.at(cluster.size()-1);
  auto it = abc.begin();
  int curr_time = (*it)->GetTimeSlice();
  int time = cell.GetTimeSlice();

//   int curr_face = (*it)->get_allcell().at(0)->get_face();
//   int face = cell.get_allcell().at(0)->get_face();
//   if (curr_face != face) return 0;
//   //std::cout << time << " " << curr_time << std::endl;

  if (time - curr_time > 2){ // there is an empty time slot
    return 0; // do not add this cell
  }else if (time-curr_time==2 && offset==2){
    for (auto it2 = abc.begin(); it2!=abc.end(); ++ it2){
      SlimMergeGeomCell *mcell = (*it2);
      if (mcell->Overlap_fast(&cell,1)){
 	//add this new cell
 	std::set<SlimMergeGeomCell*>  abc2;
 	abc2.insert(&cell);
 	cluster.push_back(abc2);
	return 1;
      }
    }
    return 0;
  }else if (time - curr_time == 0 ){ // same time slice
    // need to check previous slice to judge whether to add this cell
    if (cluster.size() < 2){
      return 0;
    }else{
      std::set<SlimMergeGeomCell*>&  abc1 = cluster.at(cluster.size()-2);
      auto it1 = abc1.begin();
      int curr_time1 = (*it1)->GetTimeSlice();
      if (time - curr_time1 == 1){
 	//start to judge whether to add this cell
 	for (auto it2 = abc1.begin(); it2!=abc1.end(); ++ it2){
 	  SlimMergeGeomCell *mcell = (*it2);
	  SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)&cell;
	  if (mcell->Overlap_fast(mcell1,offset)){
	    //add this new cell
 	    cluster.at(cluster.size()-1).insert(mcell1);
	    
	    return 1;
 	  }
 	}
	return 0;
      }else if (time-curr_time1 == 2 && offset==2){
	//start to judge whether to add this cell
 	for (auto it2 = abc1.begin(); it2!=abc1.end(); ++ it2){
 	  SlimMergeGeomCell *mcell = (*it2);
	  SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)&cell;
	  if (mcell->Overlap_fast(mcell1,1)){
	    //add this new cell
 	    cluster.at(cluster.size()-1).insert(mcell1);
	    return 1;
 	  }
 	}
	return 0;
      }else{
	return 0;
      }
    }
  }else if (time - curr_time == 1 ){
    // need to start to judge whether to add this cell
    for (auto it2 = abc.begin(); it2!=abc.end(); ++ it2){
      SlimMergeGeomCell *mcell = (*it2);
      if (mcell->Overlap_fast(&cell,offset)){
 	//add this new cell
 	std::set<SlimMergeGeomCell*>  abc2;
 	abc2.insert(&cell);
 	cluster.push_back(abc2);
	return 1;
      }
    }
    return 0;
  }else{
    return 0; // something wrong !!
  }
}


void Slim3DCluster::Form_maps(int offset, std::map<const GeomCell*, GeomCellSelection>& front_cell_map, std::map<const GeomCell*, GeomCellSelection>& back_cell_map){

  //  front_cell_map.clear();
  //  back_cell_map.clear();
  
  for (size_t i=0; i<int(cluster.size())-1; i++){
    for (auto it1 = cluster.at(i).begin(); it1!= cluster.at(i).end();it1++){
      SlimMergeGeomCell *mcell1 = (*it1);
      for (auto it2 = cluster.at(i+1).begin(); it2!= cluster.at(i+1).end();it2++){
	SlimMergeGeomCell *mcell2 = (*it2);
	if (mcell1->Overlap_fast(mcell2,offset)){
	  // add things to the map
	  if (front_cell_map.find(mcell1)==front_cell_map.end()){
	    GeomCellSelection cells;
	    cells.push_back(mcell2);
	    front_cell_map[mcell1] = cells;
	  }else{
	    front_cell_map[mcell1].push_back(mcell2);
	  }
	  
	  if (back_cell_map.find(mcell2)==back_cell_map.end()){
	    GeomCellSelection cells;
	    cells.push_back(mcell1);
	    back_cell_map[mcell2] = cells;
	  }else{
	    back_cell_map[mcell2].push_back(mcell1);
	  }
	  
	}
      }
    }
  }

  //std::cout << front_cell_map.size() << " " << back_cell_map.size() << std::endl;
}


void Slim3DCluster::MergeCluster(Slim3DCluster& cluster_to_merge){
  //need to write a merge cluster alg. 
  SlimMergeCellCluster& cluster1 = cluster_to_merge.get_ordercell();
  SlimMergeCellCluster cluster2;
  
  int i =0,  j=0;
  
  int flag = 1;
  while(flag){
    {
      std::set<SlimMergeGeomCell*>&  abc = cluster.at(i);
      auto it = abc.begin();
      int time = (*it)->GetTimeSlice();
      
      std::set<SlimMergeGeomCell*>&  abc1 = cluster1.at(j);
      auto it1 = abc1.begin();
      int time1 = (*it1)->GetTimeSlice();
      
      
      if (time == time1){
	std::set<SlimMergeGeomCell*> abc2;
	//merge both and push
	abc2.insert(abc.begin(),abc.end());
	abc2.insert(abc1.begin(),abc1.end());
	cluster2.push_back(abc2);
	i++;
	j++;
      }else if (time < time1){
	// push time
	std::set<SlimMergeGeomCell*> abc2;
	abc2.insert(abc.begin(),abc.end());
	cluster2.push_back(abc2);
	i++;
      }else if (time > time1){
	// push time1
	std::set<SlimMergeGeomCell*> abc2;
	abc2.insert(abc1.begin(),abc1.end());
	cluster2.push_back(abc2);
	j++;
      }
    }
    
    //judge i and j
    if (i==cluster.size() && j==cluster1.size()){
      flag = 0;
    }else if (i==cluster.size() && j!=cluster1.size()){
      for (int k=j;k!=cluster1.size();k++){
	std::set<SlimMergeGeomCell*> abc2;
	std::set<SlimMergeGeomCell*>&  abc1 = cluster1.at(k);
	abc2.insert(abc1.begin(),abc1.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }else if (i!=cluster.size() && j==cluster1.size()){
      for (int k=i;k!=cluster.size();k++){
	std::set<SlimMergeGeomCell*> abc2;
	std::set<SlimMergeGeomCell*>&  abc = cluster.at(k);
	abc2.insert(abc.begin(),abc.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }
  }
  
  
  cluster = cluster2;
}




