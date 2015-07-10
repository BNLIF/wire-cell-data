#ifndef GeomWireCellData_Wire_h
#define GeomWireCellData_Wire_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <set>
#include <vector>
#include <iostream>

namespace WireCell {


    /** WireCell::GeomWire - geometry information about one wire.

	Any detector application of wire cell must provide wire
	information for this class and satisfy the requirements
	given in the comments on its elements.
    */
    class GeomWire {

    public:
	GeomWire(int ident = -1,
	     WirePlaneType_t plane = kUnknownWirePlaneType,  
	     int index = -1,
	     int channel = -1,
	     const Point& point1 = Point(),
	     const Point& point2 = Point());
	~GeomWire();

	/// Detector-dependent, globally unique ID number.  Negative
	/// is illegal, not guaranteed consecutive.
	int ident() const { return _ident; }
	/// The plane/direction enum of the wire 
	WirePlaneType_t plane() const { return _plane; }
	int iplane() const { return static_cast<int>(_plane); }
	/// Consecutive, zero-based index into ordered sequence of wires in their plane
	int index() const { return _index; }
	/// Detector-dependent electronics channel number, negative is illegal.
	int channel() const { return _channel; }

	/// Return first end point the wire, in System Of Units
	Point point1() const { return _point1; }
	/// Return second end point the wire, in System Of Units
	Point point2() const { return _point2; }

	/// Return the plane+index pair.
	WirePlaneIndex plane_index() const { return WirePlaneIndex(_plane, _index); }

	bool operator<(const GeomWire& rhs) const {
	    if (_plane < rhs._plane) return true;
	    return _index < rhs._index;
	}


    protected:
	int _ident;
	WirePlaneType_t _plane;
	int _index;
	int _channel;
	Point _point1, _point2;

        friend std::ostream & operator<<(std::ostream &os, const GeomWire& gw);
    };

    std::ostream & operator<<(std::ostream &os, const GeomWire& gw);

    typedef std::pair<const GeomWire*, const GeomWire*> GeomWirePair;

    /// Used to store definitive, ordered set of wires.
    //  note: operator<() is defined in the class.
    typedef std::set<GeomWire> GeomWireSet;
    
    /// Used to temporarily construct some non-owning sub-set of cells
    typedef std::vector<const GeomWire*> GeomWireSelection;

    typedef std::map<const GeomWire*, float> WireChargeMap; 

    typedef std::map<const GeomWire*, int> WireIndexMap;

    /// Sort a GeomWireSelection by plane+index
    void sort_by_planeindex(GeomWireSelection& ws);

    /// Sort a GeomWireSelection by channel
    void sort_by_channel(GeomWireSelection& ws);
	
} // namespace WireCell
#endif
