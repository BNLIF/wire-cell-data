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
  
  if (u_proj!=0) delete u_proj;
  if (v_proj!=0) delete v_proj;
  if (w_proj!=0) delete w_proj;
}

GeomCellSelection Slim3DCluster::get_allcell(){
  gcluster.clear();
  for (int i=0;i!=cluster.size();i++){
    std::set<SlimMergeGeomCell*>  abc = cluster.at(i);
    for (auto it = abc.begin(); it!=abc.end(); ++ it){
      gcluster.push_back(*it);
    }
  }
  return gcluster;
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


int Slim3DCluster::AddCell(SlimMergeGeomCell &cell,int offset){
  std::set<SlimMergeGeomCell*>  abc = cluster.at(cluster.size()-1);
  auto it = abc.begin();
  int curr_time = (*it)->GetTimeSlice();
  int time = cell.GetTimeSlice();

//   int curr_face = (*it)->get_allcell().at(0)->get_face();
//   int face = cell.get_allcell().at(0)->get_face();
//   if (curr_face != face) return 0;
//   //std::cout << time << " " << curr_time << std::endl;


  if (time - curr_time > 1){ // there is an empty time slot
    return 0; // do not add this cell
  }else if (time - curr_time == 0 ){ // same time slice
    // need to check previous slice to judge whether to add this cell
    if (cluster.size() < 2){
      return 0;
    }else{
      std::set<SlimMergeGeomCell*>  abc1 = cluster.at(cluster.size()-2);
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



void Slim3DCluster::MergeCluster(Slim3DCluster& cluster_to_merge){
  //need to write a merge cluster alg. 
  SlimMergeCellCluster cluster1 = cluster_to_merge.get_ordercell();
  SlimMergeCellCluster cluster2;
  
  int i =0,  j=0;
  
  int flag = 1;
  while(flag){
    std::set<SlimMergeGeomCell*>  abc = cluster.at(i);
    auto it = abc.begin();
    int time = (*it)->GetTimeSlice();
      
    std::set<SlimMergeGeomCell*>  abc1 = cluster1.at(j);
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

    //judge i and j
    if (i==cluster.size() && j==cluster1.size()){
      flag = 0;
    }else if (i==cluster.size() && j!=cluster1.size()){
      for (int k=j;k!=cluster1.size();k++){
	std::set<SlimMergeGeomCell*> abc2;
	std::set<SlimMergeGeomCell*>  abc1 = cluster1.at(k);
	abc2.insert(abc1.begin(),abc1.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }else if (i!=cluster.size() && j==cluster1.size()){
      for (int k=i;k!=cluster.size();k++){
	std::set<SlimMergeGeomCell*> abc2;
	std::set<SlimMergeGeomCell*>  abc = cluster.at(k);
	abc2.insert(abc.begin(),abc.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }
  }
  
  
  cluster = cluster2;
}




