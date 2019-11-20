#include "WCPData/SpaceCell.h"
#include <vector>
using namespace WCP;

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



SpaceCell::SpaceCell(int ncluster, const GeomCell& cell1, double x, double q, double thickness)
  : cell(&cell1),
    _q(q),
    _x(x),
    _thickness(thickness)
    // _ncluster(ncluster)
{
}
