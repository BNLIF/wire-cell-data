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
  return 0;
}

