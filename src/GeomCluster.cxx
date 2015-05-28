#include "WireCellData/GeomCluster.h"

using namespace std;
using namespace WireCell;

void GeomCluster::MergeCluster(GeomCluster& cluster_to_merge){
  //need to write a merge cluster alg. 
  MergeCellCluster cluster1 = cluster_to_merge.get_ordercell();
  MergeCellCluster cluster2;
  
  int i =0,  j=0;
  
  int flag = 1;
  while(flag){
    std::set<const MergeGeomCell*>  abc = cluster.at(i);
    auto it = abc.begin();
    int time = (*it)->GetTimeSlice();
      
    std::set<const MergeGeomCell*>  abc1 = cluster1.at(j);
    auto it1 = abc1.begin();
    int time1 = (*it1)->GetTimeSlice();
   
    
    if (time == time1){
      std::set<const MergeGeomCell*> abc2;
      //merge both and push
      abc2.insert(abc.begin(),abc.end());
      abc2.insert(abc1.begin(),abc1.end());
      cluster2.push_back(abc2);
      i++;
      j++;
    }else if (time < time1){
      // push time
      std::set<const MergeGeomCell*> abc2;
      abc2.insert(abc.begin(),abc.end());
      cluster2.push_back(abc2);
      i++;
    }else if (time > time1){
      // push time1
      std::set<const MergeGeomCell*> abc2;
      abc2.insert(abc1.begin(),abc1.end());
      cluster2.push_back(abc2);
      j++;
    }

    //judge i and j
    if (i==cluster.size() && j==cluster1.size()){
      flag = 0;
    }else if (i==cluster.size() && j!=cluster1.size()){
      for (int k=j;k!=cluster1.size();k++){
	std::set<const MergeGeomCell*> abc2;
	std::set<const MergeGeomCell*>  abc1 = cluster1.at(k);
	abc2.insert(abc1.begin(),abc1.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }else if (i!=cluster.size() && j==cluster1.size()){
      for (int k=i;k!=cluster.size();k++){
	std::set<const MergeGeomCell*> abc2;
	std::set<const MergeGeomCell*>  abc = cluster.at(k);
	abc2.insert(abc.begin(),abc.end());
	cluster2.push_back(abc2);
      }
      flag = 0;
    }
  }
  
  
  cluster = cluster2;
}

GeomCellSelection GeomCluster::get_allcell(){
  gcluster.clear();
  for (int i=0;i!=cluster.size();i++){
    std::set<const MergeGeomCell*>  abc = cluster.at(i);
    for (auto it = abc.begin(); it!=abc.end(); ++ it){
      gcluster.push_back(*it);
    }
  }
  return gcluster;
}

GeomCluster::GeomCluster(MergeGeomCell &cell){
  std::set<const MergeGeomCell*>  abc;
  abc.insert(&cell);
  cluster.push_back(abc);
}

GeomCluster::~GeomCluster(){
}

int GeomCluster::AddCell(MergeGeomCell &cell){
  std::set<const MergeGeomCell*>  abc = cluster.at(cluster.size()-1);
  auto it = abc.begin();
  int curr_time = (*it)->GetTimeSlice();
  int time = cell.GetTimeSlice();

  if (time - curr_time > 1){ // there is an empty time slot
    return 0; // do not add this cell
  }else if (time - curr_time == 0 ){ // same time slice
    // need to check previous slice to judge whether to add this cell
    if (cluster.size() < 2){
      return 0;
    }else{
      std::set<const MergeGeomCell*>  abc1 = cluster.at(cluster.size()-2);
      auto it1 = abc1.begin();
      int curr_time1 = (*it1)->GetTimeSlice();
      if (time - curr_time1 == 1){
	//start to judge whether to add this cell
	for (auto it2 = abc1.begin(); it2!=abc1.end(); ++ it2){
	  const MergeGeomCell *mcell = (*it2);
	  if (mcell->Overlap(cell)){
	    //add this new cell
	    cluster.at(cluster.size()-1).insert(&cell);
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
      const MergeGeomCell *mcell = (*it2);
      if (mcell->Overlap(cell)){
	//add this new cell
	std::set<const MergeGeomCell*>  abc2;
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
