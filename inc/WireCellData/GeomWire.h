#ifndef GeomWireCellData_Wire_h
#define GeomWireCellData_Wire_h

#include "WireCellData/Units.h"
#include "WireCellData/Point.h"

#include <set>
#include <vector>
#include <iostream>

namespace WireCell {

    /// Wire set plane/direction types
  enum WirePlaneType_t {kUwire, kVwire, kYwire, kU1wire, kV1wire, kY1wire, kUnknownWirePlaneType = -1};

    /// A pair of wire plane/direction type and index w/in that plane of wires
    typedef std::pair<WirePlaneType_t, int> WirePlaneIndex;


    /** WireCell::GeomWire - geometry information about one wire.

	Any detector application of wire cell must provide wire
	information for this class and satisfy the requirements
	given in the comments on its elements.
    */
    class GeomWire {

    public:
      /*
        GeomWire(int ident = -1,
		 WirePlaneType_t plane = kUnknownWirePlaneType,  
		 int index = -1,
		 int channel = -1,
		 const Point& point1 = Point(),
		 const Point& point2 = Point(),
		 char segment = 0,
		 char face = 0,
		 short apa = 0);
      */
       GeomWire(unsigned int ident = 0,
		WirePlaneType_t plane = kUnknownWirePlaneType,  
		int index = -1,
		int channel = -1,
		const Point& point1 = Point(),
		const Point& point2 = Point(),
		char segment = 0,
		int face = -999,
		short apa = 0,
		short cryo = -1);
       
       ~GeomWire();

	/// Detector-dependent, globally unique ID number.  Negative
	/// is illegal, not guaranteed consecutive.
	int ident() const { return _ident; }
	/// The plane/direction enum of the wire 
	WirePlaneType_t plane() const { return _plane; }
	int iplane() const { return (int)_plane; }

	/// Consecutive, zero-based index into ordered sequence of wires in their plane
	int index() const { return _index; }
	/// Detector-dependent electronics channel number, negative is illegal.
	int channel() const { return _channel; }

	/// Return first end point the wire, in System Of Units
	Point point1() const { return _point1; }
	/// Return second end point the wire, in System Of Units
	Point point2() const { return _point2; }

	/// Return the number of wire segments between the channel
	/// input and this wire.  Wire directly attached to channel
	/// input is segment==0.
	int segment() const { return _segment; }

	/// Return which side of the APA this wire exists.  Nominal or
	/// "front" face is 0, other or "back" face is 1.
	int face() const { return _face; }

	/// Return the APA number associated with this wire.
	int apa() const { return _apa; }

	// Return which cryostat this wire resides.
	int cryo() const { return _cryo; }

	/// Return the plane+index pair.
	WirePlaneIndex plane_index() const { return WirePlaneIndex(_plane, _index); }

	bool operator<(const GeomWire& rhs) const {
	  if (cryo() == rhs.cryo()) {
	    if (apa() == rhs.apa()) {
	      if (face() == rhs.face()) {
		if (plane() == rhs.plane()) {
		  return index() < rhs.index();
		}
		return plane() < rhs.plane();
	      }
	      return face() < rhs.face();
	    }
	    return apa() < rhs.apa();
	  }
	  return cryo() < rhs.cryo();
	}

    protected:
	unsigned int _ident;
	WirePlaneType_t _plane;
	int _index;
	int _channel;
	Point _point1, _point2;
	char _segment;
	int _face;
	short _apa;
	short _cryo;

        friend std::ostream & operator<<(std::ostream &os, const GeomWire& gw);
    };

    std::ostream & operator<<(std::ostream &os, const GeomWire& gw);


    struct GeomWireComparePtr {
	bool operator() (const GeomWire* a, const GeomWire* b) const {
	    if (a && b) {
		return (*a) < (*b);
	    }
	    return false;
	}
    };

    typedef std::pair<const GeomWire*, const GeomWire*> GeomWirePair;

    /// Used to store definitive, ordered set of wires
    typedef std::set<GeomWire> GeomWireSet;
    
    typedef std::set<const GeomWire*, GeomWireComparePtr> GeomWirePtrSet;

    /// Used to temporarily construct some sub-set of cells
    typedef std::vector<const GeomWire*> GeomWireSelection;

    typedef std::map<const GeomWire*, float> WireChargeMap; 

    typedef std::map<const GeomWire*, int> WireIndexMap;

    /// Sort a GeomWireSelection by plane+index
    void sort_by_planeindex(GeomWireSelection& ws);

    /// Sort a GeomWireSelection by wire id
    void sort_by_ident(GeomWireSelection& ws);
    
    /// Sort a GeomWireSelection by channel
    void sort_by_channel(GeomWireSelection& ws);
	
    typedef std::vector<int> WireSelection;
    typedef std::vector<WireSelection> WireSelectionV;
    typedef std::map<int,WireSelection> WireMap;

    typedef std::map<int,std::pair<int,int>> ChirpMap;

} // namespace WireCell
#endif
