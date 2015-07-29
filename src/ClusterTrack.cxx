#include "WireCellData/ClusterTrack.h"
#include "TVector3.h"

using namespace WireCell;


ClusterTrack::ClusterTrack(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
  hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
}

ClusterTrack::~ClusterTrack(){
  delete hough;
}


void ClusterTrack::AddMSCell(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
}


void ClusterTrack::SC_Hough(Point& p, float dis){
  hough->Reset();
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
      q = cell->q();
      
      
      
      TVector3 vec(x-x0,y-y0,z-z0);
      // sc_theta.push_back(vec.Theta());
      // sc_phi.push_back(vec.Phi());
      // sc_q.push_back(q);
      if (dis <= 0){
	hough->Fill(vec.Theta(),vec.Phi(),q);
      }else{
	if (sqrt(pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2))<dis)
	  hough->Fill(vec.Theta(),vec.Phi(),q);
      }
    }
  }
}

float ClusterTrack::Get_Theta(){
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  return hough->GetXaxis()->GetBinCenter(a);
}


float ClusterTrack::Get_Phi(){
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  return hough->GetYaxis()->GetBinCenter(b);
}
