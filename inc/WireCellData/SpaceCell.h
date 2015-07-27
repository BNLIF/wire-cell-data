#ifndef SpaceWireCellData_Cell_h
#define SpaceWireCellData_Cell_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "WireCellData/GeomCell.h"

#include <set>
#include <vector>
#include <iostream>
#include <list>

namespace WireCell {

    /** WireCell::Cell - information about one space cell
     */
    class SpaceCell {
    public:
      SpaceCell(int ncluster, const GeomCell& cell, double x, double q, double thickness);

      ~SpaceCell();
	
      double x(){return _x;};
      double y(){return _y;};
      double z(){return _z;};
      double q(){return _q;};
      double area(){return _area;};
      int ncluster(){return _ncluster;};
      double thickness(){return _thickness;};
      PointVector& boundary(){return _boundary;};
	
    protected:
      double _x,_y,_z,_q;
      double _thickness;
      double _area;
      PointVector _boundary;
      int _ncluster;
      
};
    
    /// Used to temporarily collect some subset
    typedef std::vector<const WireCell::SpaceCell*> SpaceCellSelection;

}
#endif
