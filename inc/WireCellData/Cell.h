#ifndef WireCellData_Cell_h
#define WireCellData_Cell_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <list>
#include <vector>

namespace WireCellData {

    /** WireCellData::Cell - information about one cell
     */
    struct Cell {

	Cell(int ident = 0, float area=0.0 * units::centimeter2, 
	     const Point& center=Point(), 
	     const PointVector& boundary = PointVector());

	~Cell();
	
	/// A globally unique ID number (0 is illegal)
	int ident;
	/// The cross sectional area in [length^2] of the cell
	float area;
	/// A representative center point  
	Point center;
	/// A list of point giving the outline of the cell.
	PointVector boundary;
    };

    /// Used to store definitive set of cells
    typedef std::list<WireCellData::Cell> CellSet;
    /// Used to record some view into the set of cells
    typedef std::vector<WireCellData::Cell*> CellSelection;

}
#endif
