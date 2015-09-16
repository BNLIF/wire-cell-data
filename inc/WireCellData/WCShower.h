#ifndef WCShower_h
#define WCShower_h

#include "WireCellData/WCVertex.h"
#include "WireCellData/MergeSpaceCell.h"

namespace WireCell{
  
  class WCShower{
  public:
    WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cells,MergeSpaceCellMap& mcells_map);
    ~WCShower();
	    
    MergeSpaceCellSelection& get_all_cells(){return all_cells;};

    void Iterate(MergeSpaceCell *curr_cell, WireCell::MergeSpaceCellSelection &curr_cells);
    bool IsShower();

  protected:
    WCVertex *vertex;
    WCTrack *track;
    MergeSpaceCellSelection all_cells;
    MergeSpaceCellSelection& exclude_cells;
    MergeSpaceCellMap& mcells_map;
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
  
}

#endif
