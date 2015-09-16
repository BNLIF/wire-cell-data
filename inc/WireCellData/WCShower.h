#ifndef WCShower_h
#define WCShower_h

#include "WireCellData/WCVertex.h"
#include "WireCellData/MergeSpaceCell.h"

namespace WireCell{
  
  class WCShower{
  public:
    WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cells,MergeSpaceCellMap& mcells_map);
    ~WCShower();
	    
    MergeSpaceCellSelection& get_all_cells(){return all_mcells;};

    void Iterate(MergeSpaceCell *curr_cell, WireCell::MergeSpaceCellSelection &curr_cells);
    
    void SC_Hough(Point p);

    float Get_Theta(){return theta_hough;};
    float Get_Phi(){return phi_hough;};

    bool IsShower();
    


  protected:
    WCVertex *vertex;
    WCTrack *track;
    MergeSpaceCellSelection all_mcells;
    MergeSpaceCellSelection& exclude_mcells;
    MergeSpaceCellMap& mcells_map;

    float theta_hough;
    float phi_hough;
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
  
}

#endif
