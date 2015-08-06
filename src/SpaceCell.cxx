#include "WireCellData/SpaceCell.h"
#include <vector>
using namespace WireCell;

// SpaceCell::SpaceCell(int ncluster, const GeomCell& cell, double x, double q, double thickness)
//   : _x(x),
//     _q(q),
//     _thickness(thickness),
//     _ncluster(ncluster)
// {
//   _y = cell.center().y;
//   _z = cell.center().z;
//   _area = cell.cross_section();
//   _boundary = cell.boundary();
// }



SpaceCell::SpaceCell(int ncluster, const GeomCell& cell, double x, double q, double thickness)
  : cell(&cell),
    _q(q),
    _thickness(thickness)
    // _ncluster(ncluster)
{
}
