#ifndef GeomWireCellData_Wire_h
#define GeomWireCellData_Wire_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <set>
#include <vector>
#include <iostream>

namespace WireCell {

    /// Wire set plane/direction types
    enum WirePlaneType_t {kUwire, kVwire, kYwire, kUnknownWirePlaneType = -1};

    /// A pair of wire plane/direction type and index w/in that plane of wires
    typedef std::pair<WirePlaneType_t, int> WirePlaneIndex;


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

    private:
	int _ident;
	WirePlaneType_t _plane;
	int _index;
	int _channel;
	Point _point1, _point2;

        friend std::ostream & operator<<(std::ostream &os, const GeomWire& gw);
    };

    std::ostream & operator<<(std::ostream &os, const GeomWire& gw);

    /// Compare plane+index 
    struct GeomWireCompare {
	bool operator() (const GeomWire& a, const GeomWire& b) const {
	    if (a.plane() < b.plane()) return true;
	    return a.index() < b.index();
	}
    };

    typedef std::pair<const GeomWire*, const GeomWire*> GeomWirePair;

    /// Used to store definitive, ordered set of wires
    typedef std::set<GeomWire, GeomWireCompare> GeomWireSet;

    /// Used to temporarily construct some sub-set of cells
    typedef std::vector<const GeomWire*> GeomWireSelection;

    /// Sort a GeomWireSelection by plane+index
    void sort_by_planeindex(GeomWireSelection& ws);

    /// Sort a GeomWireSelection by channel
    void sort_by_channel(GeomWireSelection& ws);
	
} // namespace WireCell
#endif
