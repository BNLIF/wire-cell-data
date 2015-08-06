#ifndef GeomWireCellData_Cell_h
#define GeomWireCellData_Cell_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"
#include "WireCellData/Edge.h"

#include <set>
#include <vector>
#include <iostream>
#include <list>

namespace WireCell {

    /** WireCell::Cell - information about one cell
     */
    class GeomCell {
    public:
	GeomCell(int ident = 0, 
		 const PointVector& boundary = PointVector());
	GeomCell(const GeomCell *cell);

	~GeomCell();
	
	/// A globally unique ID number (0 is illegal)
	int ident() const { return _ident; }
	/// The cross sectional area of cell parallel to wire plane in [length^2] 
	double cross_section() const;
	/// Unbiased "center of mass" of boundary points
	Point center() const;
	/// A list of point giving the outline of the cell.
	PointVector boundary() const { return _boundary;}
	EdgeVector edge() const {return _edge;}
	const EdgeVector* redge() const { return &_edge;}

	bool operator== (const GeomCell &b) const { return this->ident()==b.ident();}

	
    protected:
	int _ident;
	PointVector _boundary;
	EdgeVector _edge;
	int order_boundary();
	mutable int flag_center;
	mutable int flag_cross_section;
	mutable Point ret;
	mutable double area;

        friend std::ostream & operator<<(std::ostream &os, const GeomCell& gc);
    };

    std::ostream & operator<<(std::ostream &os, const GeomCell& gc);

    /// Compare ident
    struct GeomCellCompare {
      bool operator() (const GeomCell& a, const GeomCell& b) const {
	
	return a.ident() < b.ident();
      }
      
    };

    /// Used to store a definitive, ordered set of cells
    typedef std::set<WireCell::GeomCell, GeomCellCompare> GeomCellSet;
    
    /// Used to temporarily collect some subset
    typedef std::vector<const WireCell::GeomCell*> GeomCellSelection;
    typedef std::vector<WireCell::GeomCellSelection> GeomCellSelectionV;
    typedef std::list<const WireCell::GeomCell*> GeomCellList;

    typedef std::map<const GeomCell*, float> CellChargeMap; 

    typedef std::map<const GeomCell*, int> CellIndexMap;
    typedef std::map<const Edge*, const GeomCell*> EdgeCellMap;
}
#endif
