#ifndef WCShower_h
#define WCShower_h

#include "WireCellData/WCVertex.h"
#include "WireCellData/MergeSpaceCell.h"

namespace WireCell{
  
  class WCShower{
  public:
    WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cell,MergeSpaceCellMap& mcells_map);
    ~WCShower();
	    
  protected:
    WCVertex *vertex;
    WCTrack *track;
    
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
  
}

#endif
