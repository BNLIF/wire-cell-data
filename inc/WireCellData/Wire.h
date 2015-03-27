#ifndef WireCellData_Wire_h
#define WireCellData_Wire_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <list>
#include <vector>

namespace WireCell {

    /// Wire set plane/direction types
    enum WirePlaneType_t {kUtype, kVtype, kYtype, kUnknown = -1};

    /// A pair of wire plane/direction type and index w/in that plane of wires
    typedef std::pair<WirePlaneType_t, int> WirePlaneIndex;

    /** WireCell::Wire - information about one wire.

	Any detector application of wire cell must provide wire
	information with this struct and satisfy the requirements
	given in the comments on its elements.
    */
    struct Wire {

	Wire(int ident = -1,
	     WirePlaneType_t plane = kUnknown,  
	     int index = -1,
	     int channel = -1,
	     const Point& point1 = Point(),
	     const Point& point2 = Point());
	~Wire();

	/// Detector-dependent, globally unique ID number.  Negative is illegal, not guaranteed consecutive.
	int ident;
	// The plane/direction enum of the wire 
	WirePlaneType_t plane;
	// Consecutive, zero-based index into ordered sequence of wires in their plane
	int index;
	// Detector-dependent electronics channel number, negative is illegal.
	int channel;
	// End points of the wire, in System Of Units, no origin specified.
	Point point1, point2;

	WirePlaneIndex plane_index() const { return WirePlaneIndex(plane, index); }
    };

    /// Used to store definitive set of wires
    typedef std::list<WireCell::Wire> WireSet;
    /// Used to record some view into the set of cells
    typedef std::vector<const WireCell::Wire*> WireSelection;


    // In-place sorts
    void sort_by_planeindex(WireSelection& ws);
    void sort_by_channel(WireSelection& ws);
	
} // namespace WireCell
#endif
