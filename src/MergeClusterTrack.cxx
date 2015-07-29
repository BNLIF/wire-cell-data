
#include "WireCellData/MergeClusterTrack.h"

using namespace WireCell;

MergeClusterTrack::MergeClusterTrack(ClusterTrack *ctrack){
  ctracks.push_back(ctrack);
  
  // add into list ... 
  all_mcells_list.assign(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end());
}

void MergeClusterTrack::Update(){
  all_mcells.clear();
  for (auto it = all_mcells_list.begin(); it!=all_mcells_list.end(); it++){
    all_mcells.push_back(*it);
  }
}
