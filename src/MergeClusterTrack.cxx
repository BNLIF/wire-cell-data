#include "WireCellData/Singleton.h"
#include "WireCellData/TPCParams.h"

#include "WireCellData/MergeClusterTrack.h"

#include "TVector3.h"
using namespace WireCell;

MergeClusterTrack::MergeClusterTrack(MergeSpaceCellSelection& mcells){
  all_mcells_list.assign(mcells.begin(),mcells.end());
  Update();
}

MergeClusterTrack::MergeClusterTrack(ClusterTrack *ctrack){
  ctracks.push_back(ctrack);
  
  // add into list ... 
  all_mcells_list.assign(ctrack->Get_allmcells().begin(),ctrack->Get_allmcells().end());

  Update();
  
}

MergeClusterTrack::~MergeClusterTrack(){
  //delete hough;
  // for (int i=0;i!=ctracks.size();i++){
  //   delete ctracks.at(i);
  // }
  // ctracks.clear();
  
}

float MergeClusterTrack::Get_Theta(){
  return theta_hough;
}


float MergeClusterTrack::Get_Phi(){
 
  return phi_hough;
}


int MergeClusterTrack::Get_TimeLength(){
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *smcell = all_mcells.at(i);
    const MergeGeomCell *mcell = smcell->get_mcell();
    
    int time;
    if (mcell!=0){
      time = mcell->GetTimeSlice();
    }else{
      time = round(((smcell->Get_Center().x)/units::cm + 256)/(Singleton<TPCParams>::Instance().get_ts_width()/units::cm));
    }

    auto it = find(times.begin(),times.end(),time);
    if (it == times.end()){
      MergeSpaceCellSelection cells;
      cells.push_back(smcell);
      times.push_back(time);
      times_mcells.push_back(cells);
    }else{
      times_mcells.at(it-times.begin()).push_back(smcell);
    }
    
  }
  // std::set<int> times;
  // for (int i=0;i!=all_mcells.size();i++){
  //   const MergeGeomCell *mcell = all_mcells.at(i)->get_mcell();
  //   times.insert(mcell->GetTimeSlice());
  // }
  return times.size();
}

MergeSpaceCellSelection& MergeClusterTrack::Get_MSCS(int time){
  return times_mcells.at(time);
}

Point MergeClusterTrack::SC_2Hough(Point &p, float dis, int flag){
  SC_Hough(p,dis, flag);
  float theta = Get_Theta();
  float phi = Get_Phi();
  Point p1;
  p1.x = p.x + dis * sin(theta) * cos(phi);
  p1.y = p.y + dis * sin(theta) * sin(phi);
  p1.z = p.z + dis * cos(theta);

  //second hough
  SC_Hough(p1,p,dis, flag);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

  // std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << " " 
  // 	    << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;

  return p1;
}



Point MergeClusterTrack::SC_IterativeHough(Point &p, float dis, int flag){
  // do three times Hough
  //first hough
  SC_Hough(p,dis, flag);
  float theta = Get_Theta();
  float phi = Get_Phi();
  Point p1;
  p1.x = p.x + dis * sin(theta) * cos(phi);
  p1.y = p.y + dis * sin(theta) * sin(phi);
  p1.z = p.z + dis * cos(theta);

  //second hough
  SC_Hough(p1,p,dis, flag);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

  //third hough
  SC_Hough(p1,p,dis, flag);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

    
  return p1;
}

void MergeClusterTrack::SC_Hough(Point&p1, Point&p, float dis, int flag){
  TH2F* hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;

  double x,y,z,q;
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);
    if (flag == 1){
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	
	x = cell->x();
	y = cell->y();
	z = cell->z();
	//q = cell->q();
	q = 1;
	
	
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
    }else if (flag == 2){
      	x = mcell->Get_Center().x;
	y = mcell->Get_Center().y;
	z = mcell->Get_Center().z;
	//q = mcell->Get_Charge();
        q = mcell->Get_all_spacecell().size();
		
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
    }else if (flag == 3){
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	x = cell->x();
	y = cell->y();
	z = cell->z();
	//q = cell->q();
	q = 1;
	TVector3 vec(x-x1,y-y1,z-z1);
	// sc_theta.push_back(vec.Theta());
	// sc_phi.push_back(vec.Phi());
	// sc_q.push_back(q);
	if (fabs(x-x0)> Singleton<TPCParams>::Instance().get_ts_width() * 3.3){  
	  if (dis <= 0){
	    hough->Fill(vec.Theta(),vec.Phi(),q);
	  }else{
	    if (sqrt(pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2))<dis)
	      hough->Fill(vec.Theta(),vec.Phi(),q);
	  }
	}
      }
    }
  }
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta_hough =  hough->GetXaxis()->GetBinCenter(a);
  phi_hough = hough->GetYaxis()->GetBinCenter(b);

  delete hough;
}


void MergeClusterTrack::SC_Hough(Point& p, float dis, int flag){
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  //  hough->Reset();
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;
  double x,y,z,q;
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);
    if (flag == 1){
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	
	x = cell->x();
	y = cell->y();
	z = cell->z();
	//q = cell->q();
	q = 1;
	
	
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
    }else if (flag==2){
      x = mcell->Get_Center().x;
      y = mcell->Get_Center().y;
      z = mcell->Get_Center().z;
      //q = mcell->Get_Charge();
      q = mcell->Get_all_spacecell().size();
      
      
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
    }else if (flag==3){
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	
	x = cell->x();
	y = cell->y();
	z = cell->z();
	//q = cell->q();
	q = 1;
	
	
	TVector3 vec(x-x0,y-y0,z-z0);
	// sc_theta.push_back(vec.Theta());
	// sc_phi.push_back(vec.Phi());
	// sc_q.push_back(q);
	if (fabs(x-x0)>Singleton<TPCParams>::Instance().get_ts_width()*3.3){
	  if (dis <= 0){
	    hough->Fill(vec.Theta(),vec.Phi(),q);
	  }else{
	    if (sqrt(pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2))<dis)
	      hough->Fill(vec.Theta(),vec.Phi(),q);
	  }
	}
      }
    }
  }

  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta_hough =  hough->GetXaxis()->GetBinCenter(a);
  phi_hough = hough->GetYaxis()->GetBinCenter(b);
  delete hough;
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

void MergeClusterTrack::Add(MergeSpaceCell *mcell, int flag){
  // flag == -1 front
  if (flag == -1){
    all_mcells_list.push_front(mcell);
  }else if (flag == 1){
    all_mcells_list.push_back(mcell);
  }
  Update();
}

void MergeClusterTrack::AddTrack(MergeClusterTrack *mctrack){
  for (int i=0;i!=mctrack->Get_allmcells().size();i++){
    MergeSpaceCell *mcell = mctrack->Get_allmcells().at(i);
    auto it = find(all_mcells.begin(),all_mcells.end(),mcell);
    if (it == all_mcells.end()){
      all_mcells.push_back(mcell);
      all_mcells_list.push_back(mcell);
    }
  }
}

void MergeClusterTrack::Organize(){
  MergeSpaceCellSet MSC_set;
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);
    MSC_set.insert(mcell);
  }
  
  all_mcells.clear();
  all_mcells_list.clear();
  for (auto it = MSC_set.begin();it!=MSC_set.end();it++){
    all_mcells.push_back(*it);
    all_mcells_list.push_back(*it);
  }
  
}

void MergeClusterTrack::Add(MergeClusterTrack *mctrack, MergeSpaceCell *mcell1, int flag_insert_direction){
  
  int flag_loop_direction=1;
  
  auto it = find(mctrack->Get_allmcells().begin(),mctrack->Get_allmcells().end(),mcell1);
  if (it - mctrack->Get_allmcells().begin() >= mctrack->Get_allmcells().end() -1 - it){ // close to the end
    flag_loop_direction = -1;
  }
  
  for (int i=0;i!=mctrack->Get_ctracks().size();i++){
    auto it =find(ctracks.begin(),ctracks.end(),mctrack->Get_ctracks().at(i));
    if (it == ctracks.end()){
      ctracks.push_back(mctrack->Get_ctracks().at(i));
    }
  }


  

  if (flag_loop_direction == 1 && flag_insert_direction == 1){
    for (int i=0;i!=mctrack->Get_allmcells().size();i++){
      MergeSpaceCell *mcell = mctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_back(mcell);
      }
    }
  }else if (flag_loop_direction == -1 && flag_insert_direction == 1){
    for (int i=mctrack->Get_allmcells().size()-1;i>-1;i--){
      MergeSpaceCell *mcell = mctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_back(mcell);
      }
    }
  }else if (flag_loop_direction == 1 && flag_insert_direction == -1){
    for (int i=0;i!=mctrack->Get_allmcells().size();i++){
      MergeSpaceCell *mcell = mctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_front(mcell);
      }
    }
  }else{
    for (int i=mctrack->Get_allmcells().size()-1;i>-1;i--){
      MergeSpaceCell *mcell = mctrack->Get_allmcells().at(i);
      auto it = find(all_mcells_list.begin(),all_mcells_list.end(),mcell);
      if (it == all_mcells_list.end()){
	all_mcells_list.push_front(mcell);
      }
    }
  }


  Update();


  
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


bool MergeClusterTrack::CheckCell(MergeSpaceCell *scell){
  if (scell == all_mcells.front())
    return true;
  if (scell == all_mcells.back())
    return true;
  if (all_mcells.size()>=2){
    if (scell == all_mcells.at(1))
      return true;
    if (scell == all_mcells.at(all_mcells.size()-2))
      return true;
  }
  return false;
}
