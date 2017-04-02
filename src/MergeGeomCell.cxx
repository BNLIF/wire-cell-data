#include "WireCellData/MergeGeomCell.h"
#include "WireCellData/Singleton.h"
#include "WireCellData/TPCParams.h"

#include <vector>
#include <cmath>
#include <list>
#include "TVector3.h"

using namespace std;
using namespace WireCell;

void MergeGeomCell::FindCorners(GeomCellMap& cmap, GeomWireMap& wmap){
  // find edge wires
  if (flag_corner == false){
    flag_corner = true;
    GeomWirePtrSet wires_u, wires_v, wires_w;
    for (int i=0;i!=cell_all.size();i++){
      const GeomCell *cell = cell_all.at(i);
      const GeomWireSelection wires = cmap[cell];
      for (int j=0;j!=wires.size();j++){
	WirePlaneType_t plane = wires.at(j)->plane();
	if (plane==0){
	  wires_u.insert(wires.at(j));
	}else if (plane==1){
	  wires_v.insert(wires.at(j));
	}else{
	  wires_w.insert(wires.at(j));
	}
      }
    }
    edge_wires.clear();
    edge_wires.push_back(*wires_u.begin());
    auto it = wires_u.end();
    it--;
    edge_wires.push_back(*it);
    edge_wires.push_back(*wires_v.begin());
    it = wires_v.end(); 
    it --;
    edge_wires.push_back(*it);
    edge_wires.push_back(*wires_w.begin());
    it = wires_w.end();
    it --;
    edge_wires.push_back(*it);
    
    // loop through all cells, and find how many of the wires are inside the edge wires, if larger than two push, and put in counter as well.
    int wire_save[3];
    for (int i=0;i!=cell_all.size();i++){
      const GeomCell *cell = cell_all.at(i);
      const GeomWireSelection wires = cmap[cell];
      int index = 0;
      for (int j=0;j!=wires.size();j++){
	auto it = find(edge_wires.begin(),edge_wires.end(),wires.at(j));
	if (it!=edge_wires.end()){
	  wire_save[index] = it - edge_wires.begin();
	  index ++;
	}
      }
      if (index >=2){
	//std::cout << index << " " << wire_save[0] << " " << wire_save[1] << " " << wire_save[2] << std::endl;
	if (index==2){
	  corner_cells_group[wire_save[0]][wire_save[1]].push_back(cell);
	}else{
	  corner_cells_group[wire_save[0]][wire_save[1]].push_back(cell);
	  corner_cells_group[wire_save[0]][wire_save[2]].push_back(cell);
	  corner_cells_group[wire_save[1]][wire_save[2]].push_back(cell);
	}
	//	for (int j=0;j!=index;j++){
	  //
	//}
	corner_cells.push_back(cell);
	corner_cells_index[cell] = index;
	//corner_cells_index.push_back(index);
      }
    }
    
    //
  }
}

int MergeGeomCell::index(int index1, int index2){
  int val;
  if (index1 == 0 && index2 ==2){
    val = 0;
  }else if (index1 == 0 && index2 ==3){
    val = 1;
  }else if (index1 == 0 && index2 ==4){
    val = 2;
  }else if (index1 == 0 && index2 ==5){
    val = 3;
  }else if (index1 == 1 && index2 ==2){
    val = 4;
  }else if (index1 == 1 && index2 ==3){
    val = 5;
  }else if (index1 == 1 && index2 ==4){
    val = 6;
  }else if (index1 == 1 && index2 ==5){
    val = 7;
  }else if (index1 == 2 && index2 ==4){
    val = 8;
  }else if (index1 == 2 && index2 ==5){
    val = 9;
  }else if (index1 == 3 && index2 ==4){
    val = 10;
  }else if (index1 == 3 && index2 ==5){
    val = 11;
  }
  
  return val;
}
int MergeGeomCell::index1(int index){
  int val;
  if (index ==0){
    val = 0;
  }else if (index==1){
    val = 0;
  }else if (index==2){
    val = 0;
  }else if (index==3){
    val = 0;
  }else if (index==4){
    val = 1;
  }else if (index==5){
    val = 1;
  }else if (index==6){
    val = 1;
  }else if (index==7){
    val = 1;
  }else if (index==8){
    val = 2;
  }else if (index==9){
    val = 2;
  }else if (index==10){
    val = 3;
  }else if (index==11){
    val = 3;
  }
  return val;
}
int MergeGeomCell::index2(int index){
   int val;
  if (index ==0){
    val = 2;
  }else if (index==1){
    val = 3;
  }else if (index==2){
    val = 4;
  }else if (index==3){
    val = 5;
  }else if (index==4){
    val = 2;
  }else if (index==5){
    val = 3;
  }else if (index==6){
    val = 4;
  }else if (index==7){
    val = 5;
  }else if (index==8){
    val = 4;
  }else if (index==9){
    val = 5;
  }else if (index==10){
    val = 4;
  }else if (index==11){
    val = 5;
  }
  return val;
}

void MergeGeomCell::FindEdges(){
  // if (flag_edge == false){
  //   flag_edge = true;
    
    edge_wires.clear();
    
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
    
    // find the min and max wires (six wires)
    const GeomWire* uwire_min = uwires.front();
    const GeomWire* uwire_max = uwires.back();
    const GeomWire* vwire_min = vwires.front();
    const GeomWire* vwire_max = vwires.back();
    const GeomWire* wwire_min = wwires.front();
    const GeomWire* wwire_max = wwires.back();
    
    for (Int_t i=0;i!=uwires.size();i++){
      TVector3 abc1(uwires.at(i)->point1().x+uwires.at(i)->point2().x,
    		    uwires.at(i)->point1().y+uwires.at(i)->point2().y,
    		    uwires.at(i)->point1().z+uwires.at(i)->point2().z);
      TVector3 abc2(uwire_min->point1().x+uwire_min->point2().x,
    		    uwire_min->point1().y+uwire_min->point2().y,
    		    uwire_min->point1().z+uwire_min->point2().z);
      TVector3 abc3(uwire_max->point1().x+uwire_max->point2().x,
    		    uwire_max->point1().y+uwire_max->point2().y,
    		    uwire_max->point1().z+uwire_max->point2().z);
      
      if (abc1.Dot(dir_up) < abc2.Dot(dir_up)){
    	uwire_min = uwires.at(i);
      }
      if (abc1.Dot(dir_up) > abc3.Dot(dir_up)){
    	uwire_max = uwires.at(i);
      }
    }
    

    for (Int_t i=0;i!=vwires.size();i++){
      TVector3 abc1(vwires.at(i)->point1().x+vwires.at(i)->point2().x,
    		    vwires.at(i)->point1().y+vwires.at(i)->point2().y,
    		    vwires.at(i)->point1().z+vwires.at(i)->point2().z);
      TVector3 abc2(vwire_min->point1().x+vwire_min->point2().x,
    		    vwire_min->point1().y+vwire_min->point2().y,
    		    vwire_min->point1().z+vwire_min->point2().z);
      TVector3 abc3(vwire_max->point1().x+vwire_max->point2().x,
    		    vwire_max->point1().y+vwire_max->point2().y,
    		    vwire_max->point1().z+vwire_max->point2().z);
      
      if (abc1.Dot(dir_vp) < abc2.Dot(dir_vp)){
    	vwire_min = vwires.at(i);
      }
      if (abc1.Dot(dir_vp) > abc3.Dot(dir_vp)){
    	vwire_max = vwires.at(i);
      }
    }

    for (Int_t i=0;i!=wwires.size();i++){
      TVector3 abc1(wwires.at(i)->point1().x+wwires.at(i)->point2().x,
    		    wwires.at(i)->point1().y+wwires.at(i)->point2().y,
    		    wwires.at(i)->point1().z+wwires.at(i)->point2().z);
      TVector3 abc2(wwire_min->point1().x+wwire_min->point2().x,
    		    wwire_min->point1().y+wwire_min->point2().y,
    		    wwire_min->point1().z+wwire_min->point2().z);
      TVector3 abc3(wwire_max->point1().x+wwire_max->point2().x,
    		    wwire_max->point1().y+wwire_max->point2().y,
    		    wwire_max->point1().z+wwire_max->point2().z);
      
      if (abc1.Dot(dir_wp) < abc2.Dot(dir_wp)){
    	wwire_min = wwires.at(i);
      }
      if (abc1.Dot(dir_wp) > abc3.Dot(dir_wp)){
    	wwire_max = wwires.at(i);
      }
    }




    // put all the cells in these wires into edge_wires;
    edge_wires.push_back(uwire_min);
    edge_wires.push_back(uwire_max);
    edge_wires.push_back(vwire_min);
    edge_wires.push_back(vwire_max);
    edge_wires.push_back(wwire_min);
    edge_wires.push_back(wwire_max);




    // std::list<const Edge*> edgelist;
    // EdgeCellMap ecmap;
    
    // for (int i=0;i!=cell_all.size();i++){
    //   const EdgeVector* evector = cell_all[i]->redge();  
      
    //   for (int j=0;j!=evector->size();j++){
    // 	int flag = 0;
    // 	// std::cout << i << " " << j << " " << evector->at(j).first.x << " "
    // 	// 		<< evector->at(j).first.y << " " << evector->at(j).first.z << " "
    // 	// 		<< evector->at(j).second.x << " "
    // 	// 		<< evector->at(j).second.y << " " << evector->at(j).second.z << " "  
    // 	// 		<< std::endl;
    // 	for (auto it = edgelist.begin(); it!=edgelist.end();it++){
    // 	  if (CompareEdge(*(*it),evector->at(j))){
    // 	    flag = 1;
    // 	    edgelist.erase(it);
    // 	    break;
    // 	  }
    // 	}
    // 	if (flag==0){
    // 	  edgelist.push_back(&(evector->at(j)));
    // 	  ecmap[&(evector->at(j))] = cell_all[i];
    // 	} 
    //   }
    // }
    
    // for (auto it = edgelist.begin(); it!=edgelist.end();it++){
    //   auto it2 = find(edge_cells.begin(),edge_cells.end(),ecmap[*it]);
    //   if (it2==edge_cells.end()){
    // 	edge_cells.push_back(ecmap[*it]);
    //   }
    // }
    
    // if (edge_cells.size() + 25 < cell_all.size()){
    //   blob = true;
    // }else{
    //   blob = false;
    // }
    
    // std::cout << edgelist.size() << " " << ecmap.size() << " " << cell_all.size() << " " << edge_cells.size() << std::endl;
    // }
}

bool MergeGeomCell::Overlap(const MergeGeomCell &cell, float num) const{
  // FindEdges();
  // cell.FindEdges();
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


      for (int i=0;i!=cell.get_uwires().size();i++){
	//	auto it = find(uwires.begin(),uwires.end(),cell.get_uwires().at(i));
	//if (it != uwires.end()){
	//}
	for (int j=0;j!=uwires.size();j++){
	  TVector3 dir(cell.get_uwires().at(i)->point1().x - uwires.at(j)->point1().x,
		       cell.get_uwires().at(i)->point1().y - uwires.at(j)->point1().y,
		       cell.get_uwires().at(i)->point1().z - uwires.at(j)->point1().z);
	  float dis = fabs(dir.Dot(dir_up));
	  if (dis < Singleton<TPCParams>::Instance().get_pitch()*1.2){
	    flag_u = 1;
	    break;
	  }
	}
	if (flag_u == 1) break;
      }
      
      int flag_v = 0;
      for (int i=0;i!=cell.get_vwires().size();i++){
	// auto it = find(vwires.begin(),vwires.end(),cell.get_vwires().at(i));
	// if (it != vwires.end()){
	//   }
	for (int j=0;j!=vwires.size();j++){
	  TVector3 dir(cell.get_vwires().at(i)->point1().x - vwires.at(j)->point1().x,
		       cell.get_vwires().at(i)->point1().y - vwires.at(j)->point1().y,
		       cell.get_vwires().at(i)->point1().z - vwires.at(j)->point1().z);
	  float dis = fabs(dir.Dot(dir_vp));
	  if (dis < Singleton<TPCParams>::Instance().get_pitch()*1.2){
	    flag_v = 1;
	    break;
	  }
	}
	if (flag_v == 1) break;
      }

      int flag_w = 0;
      for (int i=0;i!=cell.get_wwires().size();i++){
	// auto it = find(wwires.begin(),wwires.end(),cell.get_wwires().at(i));
	// if (it != wwires.end()){
	//   }
	for (int j=0;j!=wwires.size();j++){
	  TVector3 dir(cell.get_wwires().at(i)->point1().x - wwires.at(j)->point1().x,
		       cell.get_wwires().at(i)->point1().y - wwires.at(j)->point1().y,
		       cell.get_wwires().at(i)->point1().z - wwires.at(j)->point1().z);
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
  }else{
  
    for (int i=0;i!=edge_cells.size();i++){
      const GeomCell *cell1 = edge_cells[i];
      for (int j=0;j!=cell.get_edgecells().size();j++){
	const GeomCell *cell2 = cell.get_edgecells().at(j);
	
	// for (int i=0;i!=cell_all.size();i++){
	//   const GeomCell *cell1 = cell_all[i];
	//   for (int j=0;j!=cell.get_allcell().size();j++){
	//     const GeomCell *cell2 = cell.get_allcell().at(j);
	
	Point c1 = cell1->center();
	Point c2 = cell2->center();
	//std::cout << (c1.y-c2.y)/units::cm << " " << (c1.z-c2.z)/units::cm << " " << cell_all.size() << " " << cell.get_allcell().size() << " " << i << " " << j << std::endl;
	if (fabs(c1.y-c2.y) > Singleton<TPCParams>::Instance().get_pitch()*8) continue;
	if ( fabs(c1.z-c2.z) > Singleton<TPCParams>::Instance().get_pitch()*8) continue;
	
	for (int i1=0;i1!=cell1->boundary().size();i1++){
	  Point p = (cell1->boundary())[i1];
	  for (int j1=0;j1!=cell2->boundary().size();j1++){
	    Point p1 = (cell2->boundary())[j1];
	    
	    //std::cout << p.y << " " << p.z << " " << p1.y << " " << p1.z << " " << sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m << std::endl;
	    
	    if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*num){
	      
	      return true;
	    }
	  }
	}
      }
    }

    //deal with the case where one is inside the other
    int num1 = cell_all.size();
    int num2 = cell.get_allcell().size();
    const GeomCell *cell1;
    GeomCellSelection cells;
    if (num1 < num2){
      cell1 = cell_all.at(0);
      cells = cell.get_allcell();
    }else{
      cell1 = cell.get_allcell().at(0);
      cells = cell_all;
    }
    for (int j=0;j!=cells.size();j++){
      const GeomCell *cell2 = cells.at(j);
      
      Point c1 = cell1->center();
      Point c2 = cell2->center();
      //std::cout << (c1.y-c2.y)/units::cm << " " << (c1.z-c2.z)/units::cm << " " << cell_all.size() << " " << cell.get_allcell().size() << " " << i << " " << j << std::endl;
      if (fabs(c1.y-c2.y) > Singleton<TPCParams>::Instance().get_pitch()*8) continue;
      if ( fabs(c1.z-c2.z) > Singleton<TPCParams>::Instance().get_pitch()*8) continue;
      
      for (int i1=0;i1!=cell1->boundary().size();i1++){
	Point p = (cell1->boundary())[i1];
	for (int j1=0;j1!=cell2->boundary().size();j1++){
	  Point p1 = (cell2->boundary())[j1];
	  
	  //std::cout << p.y << " " << p.z << " " << p1.y << " " << p1.z << " " << sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m << std::endl;
	  
	  if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*num){
	    
	    return true;
	  }
	}
      }
    }
    
    
    
    return false;
  }
}

int MergeGeomCell::Overlap1(const MergeGeomCell &cell, float num) const{
  int val = 0;
  for (int i=0;i!=cell_all.size();i++){
    const GeomCell *cell1 = cell_all[i];
    for (int j=0;j!=cell.get_allcell().size();j++){
      const GeomCell *cell2 = cell.get_allcell().at(j);

      Point c1 = cell1->center();
      Point c2 = cell2->center();

      if (fabs(c1.y-c2.y) > Singleton<TPCParams>::Instance().get_pitch()*3.3) continue;
      if (fabs(c1.z-c2.z) > Singleton<TPCParams>::Instance().get_pitch()*3.3) continue;
      int flag = 0;
      for (int i1=0;i1!=cell1->boundary().size();i1++){
	Point p = (cell1->boundary())[i1];
	for (int j1=0;j1!=cell2->boundary().size();j1++){
	  Point p1 = (cell2->boundary())[j1];
	  if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*num){
	    flag  = 1;
	    val ++;
	    break;
	  }
	  }
	if (flag==1) break;
      }
    }
  }
  
  return val;
}

MergeGeomCell::MergeGeomCell()
{
  flag_corner = false;
  flag_edge = false;
  blob = false;
  _ident = 0;
  
  time_slice = -1;

  contain_truth = false;

  flag_center = 0;
  flag_cross_section = 0;
  ret.x = 0;
  ret.y = 0;
  ret.z = 0;
  area = 0;
}

MergeGeomCell::MergeGeomCell(int ident, const WireCell::GeomCell& cell)
{
  flag_corner = false;
  flag_edge = false;
  blob = false;
  _ident = ident;
  //this is good, as the first cell 
  _boundary = cell.boundary();
  //_edge = cell.edge();
  //
  cell_all.push_back(&cell);
  time_slice = -1;

  uwires.clear();
  uwires.push_back(cell.get_uwire());
  vwires.clear();
  vwires.push_back(cell.get_vwire());
  wwires.clear();
  wwires.push_back(cell.get_wwire());

  contain_truth = false;
  flag_center = 0;
  flag_cross_section = 0;
  ret.x = 0;
  ret.y = 0;
  ret.z = 0;
  area = 0;
}

MergeGeomCell::MergeGeomCell(int ident, const WireCell::MergeGeomCell& cell)
{
  flag_corner = false;
  flag_edge = false;
  blob = false;
  _ident = ident;

  // this is good 
  _boundary = cell.boundary();
  //_edge = cell.edge();
  // 
  
  uwires.clear();
  uwires = cell.get_uwires();
  vwires.clear();
  vwires = cell.get_vwires();
  wwires.clear();
  wwires = cell.get_wwires();



  time_slice = cell.GetTimeSlice();

  contain_truth = cell.GetContainTruthCell();
  truth_charge = cell.GetTruthCharge();

  cell_all = cell.get_allcell();
  truth_cells = cell.get_truthcell();
}

MergeGeomCell::~MergeGeomCell(){
  cell_all.clear();
  edge_cells.clear();
  edge_wires.clear();
  corner_cells.clear();
  corner_cells_index.clear();
  truth_cells.clear();
  for (int i=0;i!=6;i++){
    for (int j=0;j!=6;j++){
      corner_cells_group[i][j].clear();
    }
  }

}

double MergeGeomCell::cross_section() const
{
  //double area = 0;
  if (flag_cross_section == 0){
    flag_cross_section = 1;
    for (int i=0;i!=cell_all.size();i++){
      area += cell_all[i]->cross_section();
    }
  }
  return area;
}

Point MergeGeomCell::center() const
{
  //Point ret(0,0,0);
  if (flag_center ==0){
    flag_center = 1;
    double sum_area = 0;
    for (int i=0;i!=cell_all.size();i++){
      Point pc = cell_all[i]->center();
      double area = fabs(cell_all[i]->cross_section());
      ret.x += pc.x * area;
      ret.y += pc.y * area;
      ret.z += pc.z * area;
      sum_area += area;
    }
    
    
    ret.x/=sum_area;
    ret.y/=sum_area;
    ret.z/=sum_area;
  }
  return ret;
}

bool MergeGeomCell::Connected(const WireCell::GeomCell& cell1,const WireCell::GeomCell& cell2){
  Point c1 = cell1.center();
  Point c2 = cell2.center();
  if (fabs(c1.y-c2.y) > Singleton<TPCParams>::Instance().get_pitch()*3.3) return false;
  if (fabs(c1.z-c2.z) > Singleton<TPCParams>::Instance().get_pitch()*3.3) return false;
  

  PointVector bd1 = cell1.boundary();
  PointVector bd2 = cell2.boundary();
  int nshare=0;
  for (int i=0;i!=bd1.size();i++){
    Point p = bd1[i];
    for (int j=0;j!=bd2.size();j++){
      Point p1 = bd2[j];
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.0003)
	nshare ++;
      if (nshare == 2)
	return true;
    }
  }
  return false;
}

void MergeGeomCell::AddNewCell(const WireCell::GeomCell& cell){
  PointVector boundary = cell.boundary();
  //EdgeVector edge = cell.edge();
  // Need to improve
  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
  //_edge.insert(_edge.end(),edge.begin(),edge.end());
  // 
  
  auto it_u = find(uwires.begin(),uwires.end(),cell.get_uwire());
  if (it_u == uwires.end()){
    uwires.push_back(cell.get_uwire());
  }

  auto it_v = find(vwires.begin(),vwires.end(),cell.get_vwire());
  if (it_v == vwires.end()){
    vwires.push_back(cell.get_vwire());
  }

  auto it_w = find(wwires.begin(),wwires.end(),cell.get_wwire());
  if (it_w == wwires.end()){
    wwires.push_back(cell.get_wwire());
  }

  cell_all.push_back(&cell);
}


int MergeGeomCell::AddCell(const WireCell::GeomCell& cell, double dis){
  // check if there are too or more shared boundary points
  // if less than two shared points, not merge
  // if there are just two, and reduce two to one, and merge
  // This part can be improved, now the boundary of the merged cell are not correct
  PointVector boundary = cell.boundary();
  //EdgeVector edge = cell.edge();
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
      
      if (fabs(p.x-p1.x) > dis) continue;
      if ( fabs(p.y-p1.y) > dis) continue; 
      if (fabs(p.z-p1.z) > dis ) continue;
      
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))<dis){
	nshare ++;
	if (nshare==2){
	  // Need to improve
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
	  //_edge.insert(_edge.end(),edge.begin(),edge.end());
	  // Need to improve
	  

	  auto it_u = find(uwires.begin(),uwires.end(),cell.get_uwire());
	  if (it_u == uwires.end()){
	    uwires.push_back(cell.get_uwire());
	  }
	  
	  auto it_v = find(vwires.begin(),vwires.end(),cell.get_vwire());
	  if (it_v == vwires.end()){
	    vwires.push_back(cell.get_vwire());
	  }
	  
	  auto it_w = find(wwires.begin(),wwires.end(),cell.get_wwire());
	  if (it_w == wwires.end()){
	    wwires.push_back(cell.get_wwire());
	  }
	  cell_all.push_back(&cell);
	  return 1;
	}
      }
      
    }
  }
 
  return 0;
    
}

int MergeGeomCell::AddCell(WireCell::MergeGeomCell& cell, double dis){

  // check if there are too or more shared boundary points
  // if less than two shared points, not merge
  // if there are just two, and reduce two to one, and merge
  // This part can be improved, now the boundary of the merged cell are not correct
  PointVector boundary = cell.boundary();
  //EdgeVector edge = cell.edge();
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
     
      if (fabs(p.x-p1.x) > dis) continue; 
      if (fabs(p.y-p1.y) > dis) continue; 
      if (fabs(p.z-p1.z) > dis) continue; 
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))<dis){
	nshare ++;
	if (nshare==2){
	  // Need to improve
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
	  //_edge.insert(_edge.end(),edge.begin(),edge.end());
	  // en

	  for (int k=0;k!=cell.get_uwires().size();k++){
	    auto it_u = find(uwires.begin(),uwires.end(),cell.get_uwires().at(k));
	    if (it_u == uwires.end()){
	      uwires.push_back(cell.get_uwires().at(k));
	    }
	  }
	  
	  for (int k=0;k!=cell.get_vwires().size();k++){
	    auto it_v = find(vwires.begin(),vwires.end(),cell.get_vwires().at(k));
	    if (it_v == vwires.end()){
	      vwires.push_back(cell.get_vwires().at(k));
	    }
	  }
	  
	  for (int k=0;k!=cell.get_wwires().size();k++){
	    auto it_w = find(wwires.begin(),wwires.end(),cell.get_wwires().at(k));
	    if (it_w == wwires.end()){
	      wwires.push_back(cell.get_wwires().at(k));
	    }
	  }


	  WireCell::GeomCellSelection temp = cell.get_allcell();
	  cell_all.insert(cell_all.end(),temp.begin(),temp.end());
	  return 1;
	}
      }
    }
  }

  return 0;
}



bool MergeGeomCell::CheckContainTruthCell(WireCell::CellChargeMap &ccmap){
  truth_charge = 0;
  for (auto it = ccmap.begin();it!=ccmap.end(); it++){
    auto itt = find(cell_all.begin(),cell_all.end(),it->first);
    if (itt!=cell_all.end()){
      truth_cells.push_back(*itt);
      truth_charge += it->second;
      contain_truth = true;
    }
  }
  return contain_truth;
}

void MergeGeomCell::Organize_edge_boundary(){
  // organize edge
  // std::list<Edge> edgelist;
  // for (int i=0;i!=_edge.size();i++){
  //   int flag = 0;
  //   for (auto it = edgelist.begin();it!=edgelist.end();it++){
  //     if (CompareEdge(*it,_edge.at(i))){
  // 	edgelist.erase(it);
  // 	flag = 1;
  // 	break;
  //     }
  //   }
  //   if (flag==0){
  //     edgelist.push_back(_edge.at(i));
  //   }
  // }
  // _edge.clear();
  // _edge.insert(_edge.begin(),edgelist.begin(),edgelist.end());

  //organize boundary
  // std::list<Point> pointlist;
  // for (int i=0;i!=_edge.size();i++){
  //   for (int j=0;j!=2;j++){
  //     Point p;
  //     if (j==0){
  // 	p = _edge.at(i).first;
  //     }else if (j==1){
  // 	p = _edge.at(i).second;
  //     }
      
  //     int flag = 0;
  //     for (auto it = pointlist.begin();it!=pointlist.end();it++){
  // 	if (ComparePoint(*it,p)){
  // 	  flag = 1;
  // 	  break;
  // 	}
  //     }
  //     if (flag==0){
  // 	pointlist.push_back(p);
  //     }
      
  //   }
  // }
  // _boundary.clear();
  // _boundary.insert(_boundary.begin(),pointlist.begin(),pointlist.end());

}

// void MergeGeomCell::Insert(const EdgeVector& edge, const PointVector& boundary){
//   // _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
//   // 	  _edge.insert(_edge.end(),edge.begin(),edge.end());
  
//   std::list<Edge> edgelist(_edge.begin(),_edge.end());
//   for (int i=0;i!=edge.size();i++){
//     int flag = 0;
//     for (auto it=edgelist.begin();it!=edgelist.end();it++){
//       if (CompareEdge(*it,edge.at(i))){
// 	edgelist.erase(it);
// 	flag = 1;
// 	break;
//       }
//     }
//     if (flag==0){
//       edgelist.push_back(edge.at(i));
//     }
//   }
//   _edge.clear();
//   _edge.insert(_edge.begin(),edgelist.begin(),edgelist.end());
  
//   // if (CompareEdge(*(*it),evector->at(j))){
//   // 	    flag = 1;
//   // 	    edgelist.erase(it);
//   // 	    break;
//   // 	  }
// }
