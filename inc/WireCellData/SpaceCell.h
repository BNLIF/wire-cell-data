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
      SpaceCell(int ncluster, GeomCell& cell, float x, float q, float thickness);

      ~SpaceCell();
	
      float x(){return _x;};
      float y(){return _y;};
      float z(){return _z;};
      float q(){return _q;};
      float area(){return _area;};
      int ncluster(){return _ncluster;};
      float thickness(){return _thickness;};
      PointVector& boundary(){return _boundary;};
	
    protected:
      float _x,_y,_z,_q;
      float _thickness;
      float _area;
      PointVector _boundary;
      int _ncluster;
      
};
    
    /// Used to temporarily collect some subset
    typedef std::vector<const WireCell::SpaceCell*> SpaceCellSelection;

}
#endif
