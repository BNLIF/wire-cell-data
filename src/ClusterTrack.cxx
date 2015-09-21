#include "WireCellData/ClusterTrack.h"
#include "TVector3.h"

using namespace WireCell;


ClusterTrack::ClusterTrack(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
  
}

ClusterTrack::~ClusterTrack(){
  // delete hough;
  // for (int i =0; i!=all_mcells.size();i++){
  //   delete all_mcells.at(i);
  // }
  // all_mcells.clear();
}

Point ClusterTrack::SC_IterativeHough(Point &p, float dis, int flag){
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
  SC_Hough(p1,p,dis,flag);
  theta = Get_Theta();
  phi = Get_Phi();
  p1.x = p1.x + dis * sin(theta) * cos(phi);
  p1.y = p1.y + dis * sin(theta) * sin(phi);
  p1.z = p1.z + dis * cos(theta);

    
  return p1;
}


int ClusterTrack::CrossNum(MergeSpaceCell *mcell1, float theta, float phi){
  int result = 0;

  float x1 = mcell1->Get_Center().x;
  float y1 = mcell1->Get_Center().y;
  float z1 = mcell1->Get_Center().z;
 
  float x2 = x1 + sin(theta) * cos(phi);
  float y2 = y1 + sin(theta) * sin(phi);
  float z2 = z1 + cos(theta);

  
 
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);

    bool result1 = false;
    
    double min_dis = 1e9;
    if (mcell == mcell1){
      result1 = true;
    }else{
      for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
	SpaceCell *cell = mcell->Get_all_spacecell().at(j);
	
	
	float x,y,z;
	x = cell->x();
	y = cell->y();
	z = cell->z();
	
	double dis;
      
	TVector3 v1(x-x1,y-y1,z-z1);
	TVector3 v2(x-x2,y-y2,z-z2);
	TVector3 v3(x2-x1,y2-y1,z2-z1);
	TVector3 v4 = v1.Cross(v2);
	dis = v4.Mag()/v3.Mag();
	
	if (dis < min_dis ) min_dis = dis;
	//std::cout << "Xin1: " << dis/units::mm << std::endl;
	
	
	if (dis < 3*units::mm){
	  result1 = true;
	  break;
	}
      }
    }
    
    //std::cout << i << " " << result1 << " " << min_dis/units::mm << std::endl;
    if (result1)
      result ++;
    
  }
  
  return result;
}




int ClusterTrack::CrossNum(Point &p, float theta, float phi){
  int result = 0;

  float x1 = p.x;
  float y1 = p.y;
  float z1 = p.z;
 
  float x2 = x1 + sin(theta) * cos(phi);
  float y2 = y1 + sin(theta) * sin(phi);
  float z2 = z1 + cos(theta);

  
 
  for (int i=0;i!=all_mcells.size();i++){
    MergeSpaceCell *mcell = all_mcells.at(i);

    bool result1 = false;
    
    double min_dis = 1e9;
  
    for (int j=0;j!=mcell->Get_all_spacecell().size();j++){
      SpaceCell *cell = mcell->Get_all_spacecell().at(j);
      

      float x,y,z;
      x = cell->x();
      y = cell->y();
      z = cell->z();

      double dis;
      
      TVector3 v1(x-x1,y-y1,z-z1);
      TVector3 v2(x-x2,y-y2,z-z2);
      TVector3 v3(x2-x1,y2-y1,z2-z1);
      TVector3 v4 = v1.Cross(v2);
      dis = v4.Mag()/v3.Mag();
      
      if (dis < min_dis ) min_dis = dis;
      //std::cout << "Xin1: " << dis/units::mm << std::endl;

     
      if (dis < 3*units::mm){
	result1 = true;
	break;
      }
    }
    
    //std::cout << i << " " << result1 << " " << min_dis/units::mm << std::endl;
    if (result1)
      result ++;
    
  }
  
  return result;
}

bool ClusterTrack::AddMSCell_anyway(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
  return true;
}

bool ClusterTrack::AddMSCell(MergeSpaceCell *cell){
  bool result = false;
  if (all_mcells.size()==1){
    all_mcells.push_back(cell);
    result = true;
  }else{
    Point p1 = all_mcells.at(all_mcells.size()-2)->Get_Center();
    Point p2 = all_mcells.at(all_mcells.size()-1)->Get_Center();
    Point p3 = cell->Get_Center();

    float theta1_old = atan((p2.y-p1.y)/(p2.x-p1.x))/3.1415926*180.;
    float theta2_old = atan((p2.z-p1.z)/(p2.x-p1.x))/3.1415926*180.;

    float theta1_new = atan((p3.y-p2.y)/(p3.x-p2.x))/3.1415926*180.;
    float theta2_new = atan((p3.z-p2.z)/(p3.x-p2.x))/3.1415926*180.;
    // if (sqrt(pow(theta1_new-theta1_old,2)+pow(theta2_new-theta1_old,2))<360.){
    //   result = true;
    //   all_mcells.push_back(cell);
    // }else{
    //   result = false;
    // }

    if (cell->Get_all_spacecell().size()<3*all_mcells.at(all_mcells.size()-1)->Get_all_spacecell().size() 
        ){
      result = true;
      all_mcells.push_back(cell);
    }else{
      result = false;
    }



  }
  return result;
}

void ClusterTrack::SC_Hough(Point&p1, Point&p, float dis, int flag){
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  
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
	// set charge into 1
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
    }else if (flag==2){
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
    }

  }
  
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta_hough = hough->GetXaxis()->GetBinCenter(a);
  phi_hough = hough->GetYaxis()->GetBinCenter(b);
  
  delete hough;
}


void ClusterTrack::SC_Hough(Point& p, float dis, int flag){
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
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
    }else if (flag == 2){
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
    }else if (flag == 3){
      if (fabs(mcell->Get_Center().x - p.x) < 0.1 * mcell->thickness()) continue;
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
    }
  }

  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta_hough = hough->GetXaxis()->GetBinCenter(a);
  phi_hough = hough->GetYaxis()->GetBinCenter(b);


  delete hough;
  
}

float ClusterTrack::Get_Theta(){
  return theta_hough;
}


float ClusterTrack::Get_Phi(){
  
  return phi_hough;
}


float ClusterTrack::Get_CTheta(Point &p){
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;

  Point p1 = all_mcells.front()->Get_Center();
  Point p2 = all_mcells.back()->Get_Center();
  
  double x,y,z;
  
  if (sqrt(pow(p1.x-x0,2)+pow(p1.y-y0,2)+pow(p1.z-z0,2))<
      sqrt(pow(p2.x-x0,2)+pow(p2.y-y0,2)+pow(p2.z-z0,2))){
    x = p2.x;
    y = p2.y;
    z = p2.z;
  }else{
    x = p1.x;
    y = p1.y;
    z = p1.z;
  }
  TVector3 vec(x-x0,y-y0,z-z0);
  return vec.Theta();
  
  
}

float ClusterTrack::Get_CPhi(Point &p){
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;

  Point p1 = all_mcells.front()->Get_Center();
  Point p2 = all_mcells.back()->Get_Center();
  
  double x,y,z;
  
  if (sqrt(pow(p1.x-x0,2)+pow(p1.y-y0,2)+pow(p1.z-z0,2))<
      sqrt(pow(p2.x-x0,2)+pow(p2.y-y0,2)+pow(p2.z-z0,2))){
    x = p2.x;
    y = p2.y;
    z = p2.z;
  }else{
    x = p1.x;
    y = p1.y;
    z = p1.z;
  }
  TVector3 vec(x-x0,y-y0,z-z0);
  return vec.Phi();
}
