#include "WireCellData/MergeGeomCell.h"

#include <vector>
#include <cmath>
using namespace std;
using namespace WireCell;

MergeGeomCell::MergeGeomCell(int ident, const WireCell::GeomCell cell)
{
  _ident = ident;
  _boundary = cell.boundary();
  cell_all.push_back(&cell);
}


int MergeGeomCell::AddCell(const WireCell::GeomCell cell){
  // check if there are too or more shared boundary points
  // if less than two shared points, not merge
  // if there are just two, and reduce two to one, and merge
  // if there are more than two, remove the duplicate points, and reduce two to one, and merge
  return 0;
}

int MergeGeomCell::AddCell(WireCell::MergeGeomCell cell){
  int flag = AddCell((WireCell::GeomCell)cell);
  if (flag==0){
    return 0;
  }else{
    WireCell::GeomCellSelection temp = cell.get_allcell();
    cell_all.insert(cell_all.end(),temp.begin(),temp.end());
    return 1;
  }
}

