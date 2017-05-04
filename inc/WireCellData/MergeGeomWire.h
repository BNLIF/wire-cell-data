#ifndef GeomWireCellData_MergeWire_h
#define GeomWireCellData_MergeWire_h

#include "WireCellData/GeomWire.h"
namespace WireCell{
  class MergeGeomWire : public WireCell::GeomWire{
  public:
    
    MergeGeomWire(int ident, const WireCell::GeomWire& wire);
    MergeGeomWire(int ident, const WireCell::MergeGeomWire& wire);
    MergeGeomWire(int ident, GeomWireSelection wires);
    MergeGeomWire(const WireCell::MergeGeomWire& wire);
    ~MergeGeomWire();

    int AddWire(const WireCell::GeomWire& wire);
    int AddWire(WireCell::MergeGeomWire& wire);
    int GetTimeSlice() const { return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};
    
    WireCell::GeomWireSelection get_allwire() const{ return wire_all;}
    
  protected:
    WireCell::GeomWireSelection wire_all;
    int time_slice;
  };

  
}

#endif
