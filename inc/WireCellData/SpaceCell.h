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

      ~SpaceCell(){};
	
      double x(){return _x;};
      double y(){return cell->center().y;};
      double z(){return cell->center().z;};
      float q(){return _q;};
      double area(){return cell->cross_section();};
      //   int ncluster(){return _ncluster;};
      float thickness(){return _thickness;};
      PointVector boundary() const{return cell->boundary();};
      const GeomCell* get_cell(){return cell;};

      /* double x(){return _x;}; */
      /* double y(){return _y;}; */
      /* double z(){return ;}; */
      /* float q(){return _q;}; */
      /* double area(){return cell->cross_section();}; */
      /* //   int ncluster(){return _ncluster;}; */
      /* float thickness(){return _thickness;}; */
      /* PointVector boundary() const{return cell->boundary();}; */
      /* const GeomCell* get_cell(){return cell;}; */
            

    protected:
      const GeomCell *cell;
      float _q;
      float _thickness;
      double _x;
      // int _ncluster;    


      /* double _x,_y,_z,_q; */
      /* double _thickness; */
      /* double _area; */
      /* PointVector _boundary; */
      /* int _ncluster; */
      
};
    
    /// Used to temporarily collect some subset
    typedef std::vector<WireCell::SpaceCell*> SpaceCellSelection;
}
#endif
