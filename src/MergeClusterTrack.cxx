
#include "WireCellData/MergeClusterTrack.h"

using namespace WireCell;

MergeClusterTrack::MergeClusterTrack(ClusterTrack *ctrack){
  ctracks.push_back(ctrack);
  
  // add into list ... 
  all_mcells_list.assign(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end());

  Update();
}

ClusterTrack* MergeClusterTrack::GetClusterTrack(MergeSpaceCell* vertex){
  ClusterTrack* result = 0;

  for (int i=0;i!=ctracks.size();i++){
    result = ctracks.at(i);
    if (result->Get_FirstMSCell() == vertex 
	|| result->Get_LastMSCell() == vertex)
      break;
  }
  
  return result;
}

void MergeClusterTrack::Add(ClusterTrack *ctrack, MergeSpaceCell *mcell1){
  ctracks.push_back(ctrack);
  int flag_insert_direction=1; 
  int flag_loop_direction=1;
  
  
  
  auto it = find(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end(),mcell1);
  if (it - ctrack->Get_allmcells().begin() >= ctrack->Get_allmcells().end() -1 - it){ // close to the end
    flag_loop_direction = -1;
  }
  
  auto it1 = find(all_mcells.begin(),all_mcells.end(),mcell1);
  if (it1 - all_mcells.begin() <= all_mcells.end() -1 - it1){
    //close to the front
    flag_insert_direction = -1;
  }

  if (flag_loop_direction == 1 && flag_insert_direction == 1){
    for (int i=0;i!=ctrack->Get_allmcells().size();i++){
      MergeSpaceCell *mcell = ctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_back(mcell);
      }
    }
  }else if (flag_loop_direction == -1 && flag_insert_direction == 1){
    for (int i=ctrack->Get_allmcells().size()-1;i>-1;i--){
      MergeSpaceCell *mcell = ctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_back(mcell);
      }
    }
  }else if (flag_loop_direction == 1 && flag_insert_direction == -1){
    for (int i=0;i!=ctrack->Get_allmcells().size();i++){
      MergeSpaceCell *mcell = ctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_front(mcell);
      }
    }
  }else{
    for (int i=ctrack->Get_allmcells().size()-1;i>-1;i--){
      MergeSpaceCell *mcell = ctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_front(mcell);
      }
    }
  }


  Update();
 
  

}


void MergeClusterTrack::Update(){
  all_mcells.clear();
  for (auto it = all_mcells_list.begin(); it!=all_mcells_list.end(); it++){
    all_mcells.push_back(*it);
  }
}
