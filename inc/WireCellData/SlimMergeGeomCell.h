#ifndef GeomWireCellData_SlimMergeCell_h
#define GeomWireCellData_SlimMergeCell_h

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWireCellMap.h"
 
namespace WireCell{
class SlimMergeGeomCell : public WireCell::GeomCell{
  public:

    SlimMergeGeomCell(int ident);
   ~SlimMergeGeomCell();

    /// Unbiased "center of mass" of boundary points
   
    void AddWire(const GeomWire *wire, WirePlaneType_t plane);
    void AddBoundary( const PointVector& boundary );

    GeomWireSelection get_uwires() const{return uwires;};
    GeomWireSelection get_vwires() const{return vwires;};
    GeomWireSelection get_wwires() const{return wwires;};

    /* int GetTimeSlice() const { return time_slice;}; */
    /* void SetTimeSlice(int time){time_slice = time;}; */

    void OrderWires();

    int GetIdent() {return _ident;};

  protected:
    int _ident;

    int order_boundary();
    WireCell::GeomWireSelection uwires;
    WireCell::GeomWireSelection vwires;
    WireCell::GeomWireSelection wwires;
     
    /* int time_slice; */
    
  };
}

#endif
