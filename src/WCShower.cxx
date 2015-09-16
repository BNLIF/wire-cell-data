#include "WireCellData/WCShower.h"

using namespace WireCell;

WCShower::WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cells,MergeSpaceCellMap& mcells_map)
  : vertex(vertex)
  , track(track)
  , exclude_cells(exclude_cells)
  , mcells_map(mcells_map)
{
  // constructor to construct things that only contain things connected to it ... 
  MergeSpaceCell *start_cell = vertex->get_msc();
  MergeSpaceCellSelection curr_cells = mcells_map[start_cell];
  Iterate(start_cell, curr_cells);
}

bool WCShower::IsShower(){
  bool result = true;
  int ncell_track = track->get_centerVP_cells().size();
  //std::cout << ncell_track << " " << all_cells.size() << std::endl;
  if (ncell_track > 0.9 * all_cells.size() && ncell_track > all_cells.size() - 8)
    return false;

  return result;
}

void WCShower::Iterate(MergeSpaceCell *curr_cell, MergeSpaceCellSelection &curr_cells){
  int flag = 0;
  auto it1 = find(all_cells.begin(),all_cells.end(),curr_cell);
    
  if (it1 == all_cells.end()){
    // Not contained, do something
    auto it3 = find(track->get_centerVP_cells().begin(), track->get_centerVP_cells().end(), curr_cell);
    if (it3 != track->get_centerVP_cells().end()){
      flag = 1;
    }else{
      auto it2 = find(exclude_cells.begin(),exclude_cells.end(),curr_cell);
      if (it2 == exclude_cells.end())
	flag = 1;
    }
  }
  
  if (flag==1){
    all_cells.push_back(curr_cell);
    for (int i=0;i!=curr_cells.size();i++){
      MergeSpaceCell *curr_cell1 = curr_cells.at(i);
      MergeSpaceCellSelection curr_cells1 = mcells_map[curr_cell1];
      Iterate(curr_cell1,curr_cells1);
    }
  }

}

WCShower::~WCShower(){
  
}
