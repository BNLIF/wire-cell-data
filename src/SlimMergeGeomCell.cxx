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

bool SlimMergeGeomCell::IsSame(SlimMergeGeomCell *mcell1){
  GeomWireSelection mcell1_uwires = mcell1->get_uwires();
  GeomWireSelection mcell1_vwires = mcell1->get_vwires();
  GeomWireSelection mcell1_wwires = mcell1->get_wwires();

  if (uwires.size() != mcell1_uwires.size() ||
      vwires.size() != mcell1_vwires.size() ||
      wwires.size() != mcell1_wwires.size()){
    return false;
  }

  for (size_t i=0;i!=uwires.size();i++){
    if (uwires.at(i)!=mcell1_uwires.at(i))
      return false;
  }
  for (size_t i=0;i!=vwires.size();i++){
    if (vwires.at(i)!=mcell1_vwires.at(i))
      return false;
  }
  for (size_t i=0;i!=wwires.size();i++){
    if (wwires.at(i)!=mcell1_wwires.at(i))
      return false;
  }
  
  
  return true;
  
}

float WireCell::SlimMergeGeomCell::Estimate_minimum_charge(){
  float u_charge = 0;
  float v_charge = 0;
  float w_charge = 0;
  
  float min_charge=1e9;

  bool flag_u, flag_v, flag_w;

  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0))==bad_planes.end()){
    for (auto it = uwires.begin(); it!= uwires.end(); it++){
      u_charge += wirechargemap[*it];
    }
    flag_u = true;
  }else{
    flag_u = false;
  }
  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1))==bad_planes.end()){
    for (auto it = vwires.begin(); it!= vwires.end(); it++){
      v_charge += wirechargemap[*it];
    }
    flag_v = true;
  }else{
    flag_v = false;
  }
  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2))==bad_planes.end()){
    for (auto it = wwires.begin(); it!= wwires.end(); it++){
      w_charge += wirechargemap[*it];
    }
    flag_w = true;
  }else{
    flag_w = false;
  }
  if (flag_u && u_charge < min_charge)     min_charge = u_charge;
  if (flag_v && v_charge < min_charge)     min_charge = v_charge;
  if (flag_w && w_charge < min_charge)     min_charge = w_charge;
  
  return min_charge;
}

float WireCell::SlimMergeGeomCell::Estimate_total_charge(){
  float total_charge = 0;
  float count = 0;
  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(0))==bad_planes.end()){
    for (auto it = uwires.begin(); it!= uwires.end(); it++){
      total_charge += wirechargemap[*it];
    }
    count ++;
  }
  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(1))==bad_planes.end()){
    for (auto it = vwires.begin(); it!= vwires.end(); it++){
      total_charge += wirechargemap[*it];
    }
    count ++;
  }
  if (find(bad_planes.begin(),bad_planes.end(),WirePlaneType_t(2))==bad_planes.end()){
    for (auto it = wwires.begin(); it!= wwires.end(); it++){
      total_charge += wirechargemap[*it];
    }
    count ++;
  }

  total_charge /=count;
  
  return total_charge;
}

float WireCell::SlimMergeGeomCell::Get_Wire_Charge(const GeomWire *wire){
  if (wirechargemap.find(wire)!=wirechargemap.end()){
    return wirechargemap[wire];
  }else{
    return 0;
  }
}

float WireCell::SlimMergeGeomCell::Get_Wire_Charge_Err(const GeomWire *wire){
  if (wirechargeerrmap.find(wire)!=wirechargeerrmap.end()){
    return wirechargeerrmap[wire];
  }else{
    return 0;
  }
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

void WireCell::SlimMergeGeomCell::AddSamplingPoints(const PointVector& sampling_points){
  sample_points = sampling_points;
}



void WireCell::SlimMergeGeomCell::AddWire(const GeomWire *wire, WirePlaneType_t plane, float charge, float charge_err){
  if (plane == WirePlaneType_t(0)){
    if (find(uwires.begin(),uwires.end(),wire)==uwires.end()){
      uwires.push_back(wire);
      wirechargemap[wire] = charge;
      wirechargeerrmap[wire] = charge_err;
    }
  }else if (plane == WirePlaneType_t(1)){
    if (find(vwires.begin(),vwires.end(),wire)==vwires.end()){
      vwires.push_back(wire);
      wirechargemap[wire] = charge;
      wirechargeerrmap[wire] = charge_err;
    }
  }else if (plane == WirePlaneType_t(2)){
    if (find(wwires.begin(),wwires.end(),wire)==wwires.end()){
      wwires.push_back(wire);
      wirechargemap[wire] = charge;
      wirechargeerrmap[wire] = charge_err;
    }
  }
}


void WireCell::SlimMergeGeomCell::OrderWires(){
  WireCell::sort_by_ident(uwires);
  WireCell::sort_by_ident(vwires);
  WireCell::sort_by_ident(wwires);
}


bool WireCell::SlimMergeGeomCell::Overlap(const WireCell::SlimMergeGeomCell* cell, float num) const{

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


bool WireCell::SlimMergeGeomCell::Overlap_fast(const WireCell::SlimMergeGeomCell* cell, int offset) const{
  // std::cout << uwires.size() << " " << vwires.size() << " " << wwires.size() << " "
  // 	    << cell->get_uwires().size() << " " << cell->get_vwires().size() << " " << cell->get_wwires().size() << std::endl;
  
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

  // std::cout << u_low_index << " " << u_high_index << " " << u1_low_index << " " << u1_high_index << " " << v_low_index << " " << v_high_index << " " << v1_low_index << " " << v1_high_index << " " << w_low_index << " " << w_high_index << " " << w1_low_index << " " << w1_high_index << " " << std::endl;
  
  if (u_low_index > u1_high_index+offset || u1_low_index > u_high_index+offset) return false;
  if (v_low_index > v1_high_index+offset || v1_low_index > v_high_index+offset) return false;
  if (w_low_index > w1_high_index+offset || w1_low_index > w_high_index+offset) return false;
  
  return true;   
}

bool WireCell::SlimMergeGeomCell::Adjacent(const WireCell::SlimMergeGeomCell* cell) const {
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

  int u_score = 0;
  int v_score = 0;
  int w_score = 0;
  
  if ( u_low_index == u1_high_index+1 || u1_low_index == u_high_index+1){
    u_score = 1;
  }else if (u_low_index <= u1_high_index && u1_low_index <= u_high_index ){
    u_score = 2;
  }
  
  if (u_score==0) return false;
  
  if ( v_low_index == v1_high_index+1 || v1_low_index == v_high_index+1){
    v_score = 1;
  }else if (v_low_index <= v1_high_index && v1_low_index <= v_high_index ){
    v_score = 2;
  }
  
  if (v_score ==0) return false;

   if ( w_low_index == w1_high_index+1 || w1_low_index == w_high_index+1){
    w_score = 1;
  }else if (w_low_index <= w1_high_index && w1_low_index <= w_high_index ){
    w_score = 2;
  }
  
  if (w_score ==0) return false;

  
  if (u_score + v_score + w_score >=5){
    return true;
  }else{
    return false;
  }
  
}


