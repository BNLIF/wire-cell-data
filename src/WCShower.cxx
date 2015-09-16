#include "WireCellData/WCShower.h"
#include "TVector3.h"

using namespace WireCell;

WCShower::WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cells,MergeSpaceCellMap& mcells_map)
  : vertex(vertex)
  , track(track)
  , exclude_mcells(exclude_cells)
  , mcells_map(mcells_map)
{
  // constructor to construct things that only contain things connected to it ... 
  MergeSpaceCell *start_cell = vertex->get_msc();
  MergeSpaceCellSelection curr_cells = mcells_map[start_cell];
  Iterate(start_cell, curr_cells);
  SC_Hough(vertex->Center());
  //  std::cout << theta_hough << " " << phi_hough << std::endl;
}


void WCShower::SC_Hough(Point p){
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  double x,y,z,q;
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);
    for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
      SpaceCell *cell = mcell->Get_all_spacecell().at(j);
      x = cell->x();
      y = cell->y();
      z = cell->z();
      q = 1;
      TVector3 vec(x-x0,y-y0,z-z0);
      hough->Fill(vec.Theta(),vec.Phi(),q);
    }
  }
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta_hough =  hough->GetXaxis()->GetBinCenter(a);
  phi_hough = hough->GetYaxis()->GetBinCenter(b);
  delete hough;
}


bool WCShower::IsShower(){
  bool result = true;
  int ncell_track = track->get_centerVP_cells().size();
  //std::cout << ncell_track << " " << all_mcells.size() << std::endl;
  if (ncell_track > 0.9 * all_mcells.size() && ncell_track > all_mcells.size() - 8)
    return false;

  return result;
}

void WCShower::Iterate(MergeSpaceCell *curr_cell, MergeSpaceCellSelection &curr_cells){
  int flag = 0;
  auto it1 = find(all_mcells.begin(),all_mcells.end(),curr_cell);
    
  if (it1 == all_mcells.end()){
    // Not contained, do something
    auto it3 = find(track->get_centerVP_cells().begin(), track->get_centerVP_cells().end(), curr_cell);
    if (it3 != track->get_centerVP_cells().end()){
      flag = 1;
    }else{
      auto it2 = find(exclude_mcells.begin(),exclude_mcells.end(),curr_cell);
      if (it2 == exclude_mcells.end())
	flag = 1;
    }
  }
  
  if (flag==1){
    all_mcells.push_back(curr_cell);
    for (int i=0;i!=curr_cells.size();i++){
      MergeSpaceCell *curr_cell1 = curr_cells.at(i);
      MergeSpaceCellSelection curr_cells1 = mcells_map[curr_cell1];
      Iterate(curr_cell1,curr_cells1);
    }
  }

}

WCShower::~WCShower(){
  
}
