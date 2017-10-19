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
   
   void AddWire(const GeomWire *wire, WirePlaneType_t plane, float charge=0, float charge_err = 0);
   void AddBoundary( const PointVector& boundary );
   void AddSamplingPoints(const PointVector& sampling_points);
   //Point get_sampling_points_center();
    
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
   
   bool Adjacent(const WireCell::SlimMergeGeomCell* cell) const;
   
     
   float Get_Wire_Charge(const GeomWire *wire);
   float Get_Wire_Charge_Err(const GeomWire *wire);
   
   float Estimate_total_charge();
   float Estimate_minimum_charge();

   float get_uq(){return uq;};
   float get_vq(){return vq;};
   float get_wq(){return wq;};
   float get_udq(){return udq;};
   float get_vdq(){return vdq;};
   float get_wdq(){return wdq;};
   float get_q(){return q;}
   Point center() const;

   void set_uq(float value){uq=value;};
   void set_vq(float value){vq=value;};
   void set_wq(float value){wq=value;};
   void set_udq(float value){udq=value;};
   void set_vdq(float value){vdq=value;};
   void set_wdq(float value){wdq=value;};
   void set_q(float value){q=value;};

   bool IsSame(SlimMergeGeomCell* mcell1);

   PointVector get_sampling_points(){return sample_points;};
   
  protected:
    int _ident;
    int time_slice;

    float uq, udq, vq, vdq, wq, wdq, q;

    PointVector sample_points;
    
    //int order_boundary();
    WireCell::GeomWireSelection uwires;
    WireCell::GeomWireSelection vwires;
    WireCell::GeomWireSelection wwires;
    std::vector<WirePlaneType_t> bad_planes;

    WireCell::WireChargeMap wirechargemap;
    WireCell::WireChargeMap wirechargeerrmap;
    
    
  };
 typedef std::vector<SlimMergeGeomCell*> SMGCSelection;
 typedef std::set<SlimMergeGeomCell*> SMGCSet;
 
}

#endif
