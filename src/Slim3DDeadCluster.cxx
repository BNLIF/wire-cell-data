#include "WCPData/Slim3DDeadCluster.h"

using namespace std;
using namespace WCP;

Slim3DDeadCluster::Slim3DDeadCluster(SlimMergeGeomCell &cell, int time_slice){
  GeomCellSetp cells;
  cells.insert(&cell);
  cluster[time_slice] = cells;
  gcluster.insert(&cell);
  
  std::set<int> time_slices;
  time_slices.insert(time_slice);
  mcell_time_map[&cell]=time_slices;
}

Slim3DDeadCluster::~Slim3DDeadCluster(){
  cluster.clear();
  gcluster.clear();
  mcell_time_map.clear();
}

int Slim3DDeadCluster::AddCell(SlimMergeGeomCell &cell, int time_slice, int offset){
  // ensure if time slice or time_slice - 1 or time_slice + 1 are inside the cluster
  bool flag_save = false;
  std::vector<int> time_slices;
  
  if (cluster.find(time_slice)!=cluster.end())
    time_slices.push_back(time_slice);
  if (cluster.find(time_slice-1)!=cluster.end())
    time_slices.push_back(time_slice-1);
  if (cluster.find(time_slice+1)!=cluster.end())
    time_slices.push_back(time_slice+1);

  if (time_slices.size()>0){
    for (size_t i = 0; i!= time_slices.size(); i++){
      GeomCellSetp& mcells = cluster[time_slices.at(i)];
      for (auto it = mcells.begin(); it!=mcells.end(); it++){
	SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)&cell;
	SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
	 if (mcell->Overlap_fast(mcell1,offset)){
	   flag_save = true;
	   break;
	 }
      }
      if (flag_save)
	break;
    }
    
    if (flag_save){
      gcluster.insert(&cell);
      if (cluster.find(time_slice)==cluster.end()){
	GeomCellSetp cells;
	cells.insert(&cell);
	cluster[time_slice] = cells;
      }else{
	cluster[time_slice].insert(&cell);
      }
      if (mcell_time_map.find(&cell)==mcell_time_map.end()){
	std::set<int> time_slices;
	time_slices.insert(time_slice);
	mcell_time_map[&cell]=time_slices;
      }else{
	mcell_time_map[&cell].insert(time_slice);
      }
      return 1;
    }
   
  }
  
  return 0;
}

bool Slim3DDeadCluster::Extend(int time_slice){
  if (cluster.find(time_slice-1)!=cluster.end()){

    cluster[time_slice] = cluster[time_slice-1];
    for (auto it = cluster[time_slice].begin(); it!= cluster[time_slice].end(); it++){
      SlimMergeGeomCell *mcell = (SlimMergeGeomCell*)(*it);
      mcell_time_map[mcell].insert(time_slice);
    }
    
    return true;
  }else{
    return false;
  }
}


bool Slim3DDeadCluster::IsContain(SlimMergeGeomCell &cell, int time_slice){
  if (cluster.find(time_slice)!=cluster.end() && mcell_time_map.find(&cell) != mcell_time_map.end()){
    return true;
  }else{
    return false;
  }
}


void Slim3DDeadCluster::MergeCluster(Slim3DDeadCluster &cluster1){
  std::map<int,GeomCellSetp>& cluster1_map = cluster1.get_cluster();
  for (auto it=cluster1_map.begin(); it!= cluster1_map.end(); it++){
    int time_slice = it->first;
    GeomCellSetp& mcells = it->second;
    for (auto it= mcells.begin(); it!= mcells.end(); it++){
      SlimMergeGeomCell *cell = (SlimMergeGeomCell*)(*it);
      
      gcluster.insert(cell);
      if (cluster.find(time_slice)==cluster.end()){
	GeomCellSetp cells;
	cells.insert(cell);
	cluster[time_slice] = cells;
      }else{
	cluster[time_slice].insert(cell);
      }
      if (mcell_time_map.find(cell)==mcell_time_map.end()){
	std::set<int> time_slices;
	time_slices.insert(time_slice);
	mcell_time_map[cell]=time_slices;
      }else{
	mcell_time_map[cell].insert(time_slice);
      }
    }
  }
}
