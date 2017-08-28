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
   
   void AddWire(const GeomWire *wire, WirePlaneType_t plane, float charge=0);
    void AddBoundary( const PointVector& boundary );
    
    
    GeomWireSelection get_uwires() const{return uwires;};
    GeomWireSelection get_vwires() const{return vwires;};
    GeomWireSelection get_wwires() const{return wwires;};

    int GetTimeSlice() const { return time_slice;}; 
    void SetTimeSlice(int time){time_slice = time;}; 

    void OrderWires();

    int GetIdent() {return _ident;};
    std::vector<WirePlaneType_t> get_bad_planes(){return bad_planes;};
    void add_bad_planes(WirePlaneType_t type);
    bool Overlap(const WireCell::SlimMergeGeomCell* cell, float num=0.1) const;
    bool Overlap_fast(const WireCell::SlimMergeGeomCell* cell, int offset=1) const;
    float Get_Wire_Charge(const GeomWire *wire);
    
  protected:
    int _ident;

    //int order_boundary();
    WireCell::GeomWireSelection uwires;
    WireCell::GeomWireSelection vwires;
    WireCell::GeomWireSelection wwires;
    std::vector<WirePlaneType_t> bad_planes;

    WireCell::WireChargeMap wirechargemap;
    int time_slice; 
    
  };
}

#endif
