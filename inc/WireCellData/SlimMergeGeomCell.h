#ifndef GeomWireCellData_SlimMergeCell_h
#define GeomWireCellData_SlimMergeCell_h

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWireCellMap.h"
 
namespace WireCell{
  class SlimMergeGeomCell : public WireCell::GeomCell{
  public:

    SlimMergeGeomCell();
    ~SlimMergeGeomCell();


    void AddWire(const GeomWire *wire, WirePlaneType_t plane);
    
    GeomWireSelection get_uwires() const{return uwires;};
    GeomWireSelection get_vwires() const{return vwires;};
    GeomWireSelection get_wwires() const{return wwires;};

    int GetTimeSlice() const { return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};

    void OrderWires();

  protected:
    WireCell::GeomWireSelection uwires;
    WireCell::GeomWireSelection vwires;
    WireCell::GeomWireSelection wwires;
     
    int time_slice;
    
  };
}

#endif
