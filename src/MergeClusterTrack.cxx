
#include "WireCellData/MergeClusterTrack.h"
#include "TVector3.h"
using namespace WireCell;

MergeClusterTrack::MergeClusterTrack(ClusterTrack *ctrack){
  ctracks.push_back(ctrack);
  
  // add into list ... 
  all_mcells_list.assign(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end());

  Update();
  hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
}

MergeClusterTrack::~MergeClusterTrack(){
  delete hough;
}

float MergeClusterTrack::Get_Theta(){
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  return hough->GetXaxis()->GetBinCenter(a);
}


float MergeClusterTrack::Get_Phi(){
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  return hough->GetYaxis()->GetBinCenter(b);
}



Point MergeClusterTrack::SC_IterativeHough(Point &p, float dis){
  // do three times Hough
  //first hough
  SC_Hough(p,dis);
  float theta = Get_Theta();
  float phi = Get_Phi();
  Point p1;
  p1.x = p.x + dis * sin(theta) * cos(phi);
  p1.y = p.y + dis * sin(theta) * sin(phi);
  p1.z = p.z + dis * cos(theta);

  //second hough
  SC_Hough(p1,p,dis);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

  //third hough
  SC_Hough(p1,p,dis);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

    
  return p1;
}

void MergeClusterTrack::SC_Hough(Point&p1, Point&p, float dis){
  hough->Reset();
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;

  double x,y,z,q;
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);
    for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
      SpaceCell *cell = mcell->Get_all_spacecell().at(j);
      
      x = cell->x();
      y = cell->y();
      z = cell->z();
      q = cell->q();
      
      
      
      TVector3 vec(x-x1,y-y1,z-z1);
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


void MergeClusterTrack::SC_Hough(Point& p, float dis){
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
