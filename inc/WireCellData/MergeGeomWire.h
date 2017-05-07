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
    void order_wires();
    
    WireCell::GeomWireSelection get_allwire() const{ return wire_all;}
    
  protected:
    bool sort_flag;
    WireCell::GeomWireSelection wire_all;
    int time_slice;
  };

  struct MergeGeomWireCompare {
    bool operator() (MergeGeomWire* a, MergeGeomWire *b) const {
      a->order_wires();
      b->order_wires();
      
      if (a->get_allwire().front()->index() < b->get_allwire().front()->index()){
	return true;
      }else if (a->get_allwire().front()->index() > b->get_allwire().front()->index()){
	return false;
      }else{
	if (a->get_allwire().back()->index() < b->get_allwire().back()->index()){
	  return true;
	}else if (a->get_allwire().back()->index() > b->get_allwire().back()->index()){
	  return false;
	}
      }
      
      return false;
    }
  };
  
  typedef std::set<MergeGeomWire*, MergeGeomWireCompare> MergeGeomWireSet;

}

#endif
