#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/TPCParams.h"
#include "TVector3.h"

using namespace WireCell;

float MergeSpaceCell::Get_Charge(){
  float sum = 0;
  for (int j=0;j!=all_spacecell.size();j++){
    SpaceCell *cell = all_spacecell.at(j);
    sum += cell->q();
  }
  return sum;
}


void MergeSpaceCell::CalMinMax(){
  max_y = -1e9;
  min_y = 1e9;
  
  max_z = -1e9;
  min_z = 1e9;
  
  for (int i=0;i!=all_spacecell.size();i++){
    SpaceCell *cell = all_spacecell.at(i);
    float y = cell->y();
    float z = cell->z();

    if (y > max_y) max_y = y;
    if (z > max_z) max_z = z;
    if (y < min_y) min_y = y;
    if (z < min_z) min_z = z;
    
  }
}


double MergeSpaceCell::ClosestDis(Point& p){
  double dis = 1e9;
  double x1 = p.x;
  double y1 = p.y;
  double z1 = p.z;
  
  for (int j=0;j!=Get_all_spacecell().size();j++){
    SpaceCell *cell = Get_all_spacecell().at(j);
    double x,y,z;
    x = cell->x();
    y = cell->y();
    z = cell->z();
    
    double dis1 = sqrt(pow(x-x1,2)+pow(y-y1,2)+pow(z-z1,2));
    if (dis1 < dis) dis = dis1;
  }

  return dis;

}

bool MergeSpaceCell::CrossCell(Point &p, float theta, float phi, int flag){
    
  float x1 = p.x;
  float y1 = p.y;
  float z1 = p.z;
 
  float x2 = x1 + sin(theta) * cos(phi);
  float y2 = y1 + sin(theta) * sin(phi);
  float z2 = z1 + cos(theta);

  
  for (int j=0;j!=Get_all_spacecell().size();j++){
    SpaceCell *cell = Get_all_spacecell().at(j);
    
    
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
    
    if (dis < Singleton<TPCParams>::Instance().get_pitch() && flag == 0){
      return true;
    }else if (dis < Singleton<TPCParams>::Instance().get_pitch() && flag == 1){
      return true;
    }
  }
  return false;
}


bool MergeSpaceCell::Overlap(MergeSpaceCell& mcell1,float num){
  if (mcell ==0){
    if (num < 0.5){
      // use the wires to determine if overlaps
      int flag_u = 0;
      
      TVector3 dir_x(1,0,0);
      TVector3 dir_u(uwires.at(0)->point1().x - uwires.at(0)->point2().x,
		     uwires.at(0)->point1().y - uwires.at(0)->point2().y,
		     uwires.at(0)->point1().z - uwires.at(0)->point2().z);
      TVector3 dir_v(vwires.at(0)->point1().x - vwires.at(0)->point2().x,
		     vwires.at(0)->point1().y - vwires.at(0)->point2().y,
		     vwires.at(0)->point1().z - vwires.at(0)->point2().z);
      TVector3 dir_w(wwires.at(0)->point1().x - wwires.at(0)->point2().x,
		     wwires.at(0)->point1().y - wwires.at(0)->point2().y,
		     wwires.at(0)->point1().z - wwires.at(0)->point2().z);

      TVector3 dir_up = dir_x.Cross(dir_u);
      TVector3 dir_vp = dir_x.Cross(dir_v);
      TVector3 dir_wp = dir_x.Cross(dir_w);

      dir_up *= 1./dir_up.Mag();
      dir_vp *= 1./dir_vp.Mag();
      dir_wp *= 1./dir_wp.Mag();


      for (int i=0;i!=mcell1.get_uwires().size();i++){
	//	auto it = find(uwires.begin(),uwires.end(),mcell1.get_uwires().at(i));
	//if (it != uwires.end()){
	//}
	for (int j=0;j!=uwires.size();j++){
	  TVector3 dir(mcell1.get_uwires().at(i)->point1().x - uwires.at(j)->point1().x,
		       mcell1.get_uwires().at(i)->point1().y - uwires.at(j)->point1().y,
		       mcell1.get_uwires().at(i)->point1().z - uwires.at(j)->point1().z);
	  float dis = fabs(dir.Dot(dir_up));
	  if (dis < Singleton<TPCParams>::Instance().get_pitch()*1.2){
	    flag_u = 1;
	    break;
	  }
	}
	if (flag_u == 1) break;
      }
      
      int flag_v = 0;
      for (int i=0;i!=mcell1.get_vwires().size();i++){
	// auto it = find(vwires.begin(),vwires.end(),mcell1.get_vwires().at(i));
	// if (it != vwires.end()){
	//   }
	for (int j=0;j!=vwires.size();j++){
	  TVector3 dir(mcell1.get_vwires().at(i)->point1().x - vwires.at(j)->point1().x,
		       mcell1.get_vwires().at(i)->point1().y - vwires.at(j)->point1().y,
		       mcell1.get_vwires().at(i)->point1().z - vwires.at(j)->point1().z);
	  float dis = fabs(dir.Dot(dir_vp));
	  if (dis < Singleton<TPCParams>::Instance().get_pitch()*1.2){
	    flag_v = 1;
	    break;
	  }
	}
	if (flag_v == 1) break;
      }

      int flag_w = 0;
      for (int i=0;i!=mcell1.get_wwires().size();i++){
	// auto it = find(wwires.begin(),wwires.end(),mcell1.get_wwires().at(i));
	// if (it != wwires.end()){
	//   }
	for (int j=0;j!=wwires.size();j++){
	  TVector3 dir(mcell1.get_wwires().at(i)->point1().x - wwires.at(j)->point1().x,
		       mcell1.get_wwires().at(i)->point1().y - wwires.at(j)->point1().y,
		       mcell1.get_wwires().at(i)->point1().z - wwires.at(j)->point1().z);
	  float dis = fabs(dir.Dot(dir_wp));
	  if (dis < Singleton<TPCParams>::Instance().get_pitch()*1.2){
	    flag_w = 1;
	    break;
	  }
	}
	if (flag_w==1) break;
      }

      // if (flag_u == 1 && flag_v == 1 && flag_w == 1){
      if (flag_u + flag_v + flag_w ==3){
	return true;
      }else{
	return false;
      }


      // // use the wires to determine if overlaps
      // int flag_u = 0;
      // for (int i=0;i!=mcell1.get_uwires().size();i++){
      // 	auto it = find(uwires.begin(),uwires.end(),mcell1.get_uwires().at(i));
      // 	if (it != uwires.end()){
      // 	  flag_u = 1;
      // 	  break;
      // 	}
      // }
      
      // int flag_v = 0;
      // for (int i=0;i!=mcell1.get_vwires().size();i++){
      // 	auto it = find(vwires.begin(),vwires.end(),mcell1.get_vwires().at(i));
      // 	if (it != vwires.end()){
      // 	  flag_v = 1;
      // 	  break;
      // 	}
      // }

      // int flag_w = 0;
      // for (int i=0;i!=mcell1.get_wwires().size();i++){
      // 	auto it = find(wwires.begin(),wwires.end(),mcell1.get_wwires().at(i));
      // 	if (it != wwires.end()){
      // 	  flag_w = 1;
      // 	  break;
      // 	}
      // }

      // if (flag_u == 1 && flag_v == 1 && flag_w == 1){
      // 	return true;
      // }else{
      // 	return false;
      // }
      
      
    }else{
      for (int i=0;i!=all_spacecell.size();i++){
	SpaceCell *cell1 = all_spacecell.at(i);
	for (int j=0;j!=mcell1.Get_all_spacecell().size();j++){
	  SpaceCell *cell2 = mcell1.Get_all_spacecell().at(j);
	  
	  
	  if (fabs(cell1->y()-cell2->y()) >  Singleton<TPCParams>::Instance().get_pitch() * 8 ) continue; 
	  if (fabs(cell1->z()-cell2->z()) > Singleton<TPCParams>::Instance().get_pitch() * 8) continue;
	  
	  for (int i1=0;i1!=cell1->boundary().size();i1++){
	    Point p = (cell1->boundary())[i1];
	    for (int j1=0;j1!=cell2->boundary().size();j1++){
	      Point p1 = (cell2->boundary())[j1];
	      if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*num){
		// std::cout << p.y << " " << p.z << " " << p1.y << " " << p1.z << std::endl;
		return true;
	      }
	    }
	  }
	}
      }
      return false;
    }
  }else{
    return mcell->Overlap(*(mcell1.get_mcell()),num);
  }
}

Point& MergeSpaceCell::Get_Center(){
  if (center_flag==1){
    return center;
  }else{
    // center.x = 0;
    // center.y = 0;
    // center.z = 0;
    
    // for (int i=0;i!=all_spacecell.size();i++){
    //   SpaceCell* cell = all_spacecell.at(i);
    //   center.x += cell->x();
    //   center.y += cell->y();
    //   center.z += cell->z();
    // }
    
    // center.x /= all_spacecell.size();
    // center.y /= all_spacecell.size();
    // center.z /= all_spacecell.size();

    double x = 0;
    double y = 0;
    double z = 0;
    for (int i=0;i!=all_spacecell.size();i++){
      SpaceCell* cell = all_spacecell.at(i);
      x += cell->x();
      y += cell->y();
      z += cell->z();
    }
    x /= all_spacecell.size();
    y /= all_spacecell.size();
    z /= all_spacecell.size();

    center.x = x;
    center.y = y;
    center.z = z;

    
    return center;
  }
}


MergeSpaceCell::~MergeSpaceCell(){
  mcell = 0;
  // for (int i=0;i!=all_spacecell.size();i++){
  //   delete all_spacecell.at(i);
  // }
  // all_spacecell.clear();
}

void MergeSpaceCell::AddSpaceCell(SpaceCell *cell){
  all_spacecell.push_back(cell);
  if (mcell == 0){
    auto it_u = find(uwires.begin(),uwires.end(),cell->get_uwire());
    if (it_u == uwires.end()){
      uwires.push_back(cell->get_uwire());
    }
    
    auto it_v = find(vwires.begin(),vwires.end(),cell->get_vwire());
    if (it_v == vwires.end()){
      vwires.push_back(cell->get_vwire());
    }
    
    auto it_w = find(wwires.begin(),wwires.end(),cell->get_wwire());
    if (it_w == wwires.end()){
      wwires.push_back(cell->get_wwire());
    }
  }
}
