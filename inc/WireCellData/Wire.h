#ifndef WireCellData_Wire_h
#define WireCellData_Wire_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <list>
#include <vector>

namespace WireCellData {

    /// Wire set plane/direction types
    enum WirePlaneType_t {kUnknown = 0, kUtype, kVtype, kYtype};

    /// A pair of wire plane/direction type and index w/in that plane of wires
    typedef std::pair<WirePlaneType_t, int> WirePlaneIndex;

    /** WireCellData::Wire - information about one wire.
    */
    struct Wire {

	Wire(int ident = 0,  	// 0 = illegal
	     WirePlaneType_t plane = kUnknown,  
	     int index = -1,	// 0-based index into one plane
	     int channel = -1,	// electronics channel
	     const Point& point1 = Point(),
	     const Point& point2 = Point());
	~Wire();

	/// globally unique ID (zero is illegal)
	int ident;
	// the plane/direction type of the wire (zero is illegal)
	WirePlaneType_t plane;
	// index into ordered sequence of wires in a plane
	int index;
	// electronics channel connected
	int channel;
	// End points of the wire
	Point point1, point2;

	WirePlaneIndex plane_index() const { return WirePlaneIndex(plane, index); }
    };

    /// Used to store definitive set of wires
    typedef std::list<WireCellData::Wire> WireSet;
    /// Used to record some view into the set of cells
    typedef std::vector<const WireCellData::Wire*> WireSelection;


    // In-place sorts
    void sort_by_planeindex(WireCellData::WireSelection& ws);
    void sort_by_channel(WireCellData::WireSelection& ws);
	
}
#endif
