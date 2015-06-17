#include "WireCellData/MergeGeomCell.h"

#include <vector>
#include <cmath>
#include <list>
using namespace std;
using namespace WireCell;

void MergeGeomCell::FindCorners(GeomCellMap& cmap, GeomWireMap& wmap){
  // find edge wires
  if (flag_corner == false){
    flag_corner = true;
    GeomWireSet1 wires_u, wires_v, wires_w;
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
	  wire_save[index] = j;
	  index ++;
	}
      }
      if (index >=2){
	for (int j=0;j!=index;j++){
	  corner_cells_group[wire_save[index]].push_back(cell);
	}
	corner_cells.push_back(cell);
	corner_cells_index[cell] = index;
	//corner_cells_index.push_back(index);
      }
    }
  }
}

void MergeGeomCell::FindEdges(){
  if (flag_edge == false){
    flag_edge = true;
    std::list<const Edge*> edgelist;
    EdgeCellMap ecmap;
    
    for (int i=0;i!=cell_all.size();i++){
      const EdgeVector* evector = cell_all[i]->redge();  
      
      for (int j=0;j!=evector->size();j++){
	int flag = 0;
	// std::cout << i << " " << j << " " << evector->at(j).first.x << " "
	// 		<< evector->at(j).first.y << " " << evector->at(j).first.z << " "
	// 		<< evector->at(j).second.x << " "
	// 		<< evector->at(j).second.y << " " << evector->at(j).second.z << " "  
	// 		<< std::endl;
	for (auto it = edgelist.begin(); it!=edgelist.end();it++){
	  if (CompareEdge(*(*it),evector->at(j))){
	    flag = 1;
	    edgelist.erase(it);
	    break;
	  }
	}
	if (flag==0){
	  edgelist.push_back(&(evector->at(j)));
	  ecmap[&(evector->at(j))] = cell_all[i];
	} 
      }
    }
    
    for (auto it = edgelist.begin(); it!=edgelist.end();it++){
      auto it2 = find(edge_cells.begin(),edge_cells.end(),ecmap[*it]);
      if (it2==edge_cells.end()){
	edge_cells.push_back(ecmap[*it]);
      }
    }
    
    if (edge_cells.size() + 9 < cell_all.size()){
      blob = true;
    }else{
      blob = false;
    }
    
    // std::cout << edgelist.size() << " " << ecmap.size() << " " << cell_all.size() << " " << edge_cells.size() << std::endl;
  }
}

bool MergeGeomCell::Overlap(const MergeGeomCell &cell) const{
  for (int i=0;i!=cell_all.size();i++){
    const GeomCell *cell1 = cell_all[i];
    for (int j=0;j!=cell.get_allcell().size();j++){
      const GeomCell *cell2 = cell.get_allcell().at(j);
      
      for (int i1=0;i1!=cell1->boundary().size();i1++){
	Point p = (cell1->boundary())[i1];
	for (int j1=0;j1!=cell2->boundary().size();j1++){
	  Point p1 = (cell2->boundary())[j1];
	  if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*4){
	    return true;
	  }
	}
      }

    }
  }
  return false;
}

MergeGeomCell::MergeGeomCell(int ident, const WireCell::GeomCell& cell)
{
  flag_corner = false;
  flag_edge = false;
  blob = false;
  _ident = ident;
  _boundary = cell.boundary();
  _edge = cell.edge();
  cell_all.push_back(&cell);
  time_slice = -1;

  contain_truth = false;
}

MergeGeomCell::MergeGeomCell(int ident, const WireCell::MergeGeomCell& cell)
{
  flag_corner = false;
  flag_edge = false;
  blob = false;
  _ident = ident;
  _boundary = cell.boundary();
  _edge = cell.edge();
  
  time_slice = cell.GetTimeSlice();

  contain_truth = cell.GetContainTruthCell();
  truth_charge = cell.GetTruthCharge();

  cell_all = cell.get_allcell();
  truth_cells = cell.get_truthcell();
}

MergeGeomCell::~MergeGeomCell(){
  cell_all.clear();
}

double MergeGeomCell::cross_section() const
{
  double area = 0;
  for (int i=0;i!=cell_all.size();i++){
    area += cell_all[i]->cross_section();
  }
  return area;
}

Point MergeGeomCell::center() const
{
  Point ret(0,0,0);
  double sum_area = 0;
  for (int i=0;i!=cell_all.size();i++){
    Point pc = cell_all[i]->center();
    double area = cell_all[i]->cross_section();
    ret.x = pc.x * area;
    ret.y = pc.z * area;
    ret.z = pc.y * area;
    sum_area += area;
  }
  
  
  ret.x/=sum_area;
  ret.y/=sum_area;
  ret.z/=sum_area;
  return ret;
}


int MergeGeomCell::AddCell(const WireCell::GeomCell& cell){
  // check if there are too or more shared boundary points
  // if less than two shared points, not merge
  // if there are just two, and reduce two to one, and merge
  // This part can be improved, now the boundary of the merged cell are not correct
  PointVector boundary = cell.boundary();
  EdgeVector edge = cell.edge();
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.0002){
  	nshare ++;
	if (nshare==2){
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
	  _edge.insert(_edge.end(),edge.begin(),edge.end());
	  cell_all.push_back(&cell);
	  return 1;
	}
      }
    }
  }
 
  return 0;
    
}

int MergeGeomCell::AddCell(WireCell::MergeGeomCell& cell){

  // check if there are too or more shared boundary points
  // if less than two shared points, not merge
  // if there are just two, and reduce two to one, and merge
  // This part can be improved, now the boundary of the merged cell are not correct
  PointVector boundary = cell.boundary();
  EdgeVector edge = cell.edge();
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.0002){
  	nshare ++;
	if (nshare==2){
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
	  _edge.insert(_edge.end(),edge.begin(),edge.end());
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
