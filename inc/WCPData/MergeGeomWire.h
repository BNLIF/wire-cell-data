#ifndef GeomWCPData_MergeWire_h
#define GeomWCPData_MergeWire_h

#include "WCPData/GeomWire.h"
namespace WCP{
  class MergeGeomWire : public WCP::GeomWire{
  public:
    
    MergeGeomWire(int ident, const WCP::GeomWire& wire);
    MergeGeomWire(int ident, const WCP::MergeGeomWire& wire);
    MergeGeomWire(int ident, GeomWireSelection wires);
    MergeGeomWire(const WCP::MergeGeomWire& wire);
    ~MergeGeomWire();

    // int ident() const { return _ident; }
    
    int AddWire(const WCP::GeomWire& wire);
    int AddWire(WCP::MergeGeomWire& wire);
    int GetTimeSlice() const { return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};
    void order_wires();
    
    WCP::GeomWireSelection get_allwire() const{ return wire_all;}
    
  protected:
    bool sort_flag;
    WCP::GeomWireSelection wire_all;
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
