#include "WireCellData/PR3DCluster.h"

using namespace WireCell;

PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
}

PR3DCluster::~PR3DCluster(){
}

// void AddCell(SlimMergeGeomCell* mcell, int *time_slices, int ntime_slice){
  
// }

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
