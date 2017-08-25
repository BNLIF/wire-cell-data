#include "WireCellData/SlimMergeGeomCell.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/TPCParams.h"

#include "TVector3.h"

using namespace WireCell;

WireCell::SlimMergeGeomCell::SlimMergeGeomCell(int ident)
  : _ident(ident)
  , time_slice(-1)
{
}



WireCell::SlimMergeGeomCell::~SlimMergeGeomCell(){
  uwires.clear();
  vwires.clear();
  wwires.clear();
}

void WireCell::SlimMergeGeomCell::add_bad_planes(WirePlaneType_t plane){
  if (find(bad_planes.begin(),bad_planes.end(),plane)!=bad_planes.end()){
  }else{
    bad_planes.push_back(plane);
  }
}

void WireCell::SlimMergeGeomCell::AddBoundary(const PointVector& boundary){
  _boundary = boundary;
  flag_center = 0;
  order_boundary();
}



void WireCell::SlimMergeGeomCell::AddWire(const GeomWire *wire, WirePlaneType_t plane){
  if (plane == WirePlaneType_t(0)){
    if (find(uwires.begin(),uwires.end(),wire)==uwires.end())
      uwires.push_back(wire);
  }else if (plane == WirePlaneType_t(1)){
    if (find(vwires.begin(),vwires.end(),wire)==vwires.end())
      vwires.push_back(wire);
  }else if (plane == WirePlaneType_t(2)){
    if (find(wwires.begin(),wwires.end(),wire)==wwires.end())
      wwires.push_back(wire);
  }
}


void WireCell::SlimMergeGeomCell::OrderWires(){
  WireCell::sort_by_ident(uwires);
  WireCell::sort_by_ident(vwires);
  WireCell::sort_by_ident(wwires);
}


bool WireCell::SlimMergeGeomCell::Overlap(const WireCell::SlimMergeGeomCell* cell, float num){

  double cut_limit = 1.2;
  
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
  

  for (int i=0;i!=cell->get_uwires().size();i++){
    //	auto it = find(uwires.begin(),uwires.end(),cell.get_uwires().at(i));
    //if (it != uwires.end()){
    //}
    for (int j=0;j!=uwires.size();j++){
      TVector3 dir(cell->get_uwires().at(i)->point1().x - uwires.at(j)->point1().x,
		   cell->get_uwires().at(i)->point1().y - uwires.at(j)->point1().y,
		   cell->get_uwires().at(i)->point1().z - uwires.at(j)->point1().z);
      float dis = fabs(dir.Dot(dir_up));
      if (dis < Singleton<TPCParams>::Instance().get_pitch()*cut_limit){
	flag_u = 1;
	break;
      }
    }
    if (flag_u == 1) break;
  }
  if (flag_u==0) return false;
      
  int flag_v = 0;
  for (int i=0;i!=cell->get_vwires().size();i++){
    // auto it = find(vwires.begin(),vwires.end(),cell.get_vwires().at(i));
    // if (it != vwires.end()){
    //   }
    for (int j=0;j!=vwires.size();j++){
      TVector3 dir(cell->get_vwires().at(i)->point1().x - vwires.at(j)->point1().x,
		   cell->get_vwires().at(i)->point1().y - vwires.at(j)->point1().y,
		   cell->get_vwires().at(i)->point1().z - vwires.at(j)->point1().z);
      float dis = fabs(dir.Dot(dir_vp));
      if (dis < Singleton<TPCParams>::Instance().get_pitch()*cut_limit){
	flag_v = 1;
	break;
      }
    }
    if (flag_v == 1) break;
  }
  if (flag_v==0) return false;  

  int flag_w = 0;
  for (int i=0;i!=cell->get_wwires().size();i++){
    // auto it = find(wwires.begin(),wwires.end(),cell.get_wwires().at(i));
    // if (it != wwires.end()){
    //   }
    for (int j=0;j!=wwires.size();j++){
      TVector3 dir(cell->get_wwires().at(i)->point1().x - wwires.at(j)->point1().x,
		   cell->get_wwires().at(i)->point1().y - wwires.at(j)->point1().y,
		   cell->get_wwires().at(i)->point1().z - wwires.at(j)->point1().z);
      float dis = fabs(dir.Dot(dir_wp));
      if (dis < Singleton<TPCParams>::Instance().get_pitch()*cut_limit){
	flag_w = 1;
	break;
      }
    }
    if (flag_w==1) break;
  }
  if (flag_w==0) return false;
  
  // if (flag_u == 1 && flag_v == 1 && flag_w == 1){
  if (flag_u + flag_v + flag_w ==3){
    return true;
  }
  
  return false;
}


bool WireCell::SlimMergeGeomCell::Overlap_fast(const WireCell::SlimMergeGeomCell* cell){
  //std::cout << uwires.size() << " " << vwires.size() << " " << wwires.size() << " "
  //	    << cell->get_uwires().size() << " " << cell->get_vwires().size() << " " << cell->get_wwires().size() << std::endl;
  
  int u_low_index = uwires.front()->index();
  int u_high_index = uwires.back()->index();

  int v_low_index = vwires.front()->index();
  int v_high_index = vwires.back()->index();

  int w_low_index = wwires.front()->index();
  int w_high_index = wwires.back()->index();

  int u1_low_index = cell->get_uwires().front()->index();
  int u1_high_index = cell->get_uwires().back()->index();

  int v1_low_index = cell->get_vwires().front()->index();
  int v1_high_index = cell->get_vwires().back()->index();

  int w1_low_index = cell->get_wwires().front()->index();
  int w1_high_index = cell->get_wwires().back()->index();

  if (u_low_index > u1_high_index || u1_low_index > u_high_index) return false;
  if (v_low_index > v1_high_index || v1_low_index > v_high_index) return false;
  if (w_low_index > w1_high_index || w1_low_index > w_high_index) return false;
  
  return true;
  
  
}
