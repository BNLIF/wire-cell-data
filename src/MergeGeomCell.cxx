#include "WireCellData/MergeGeomCell.h"

#include <vector>
#include <cmath>
using namespace std;
using namespace WireCell;



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
  _ident = ident;
  _boundary = cell.boundary();
  cell_all.push_back(&cell);
  time_slice = -1;

  contain_truth = false;
}

MergeGeomCell::MergeGeomCell(int ident, const WireCell::MergeGeomCell& cell)
{
  _ident = ident;
  _boundary = cell.boundary();
  
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
  
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.0002){
  	nshare ++;
	if (nshare==2){
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
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
  
  int nshare = 0;
  for (int i=0;i!=boundary.size();i++){
    Point p = boundary[i];
    for (int j=0;j!=_boundary.size();j++){
      Point p1 = _boundary[j];
      if (sqrt(pow(p.x-p1.x,2)+pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.0002){
  	nshare ++;
	if (nshare==2){
	  _boundary.insert(_boundary.end(),boundary.begin(),boundary.end());
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
