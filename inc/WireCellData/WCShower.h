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
    float SC_proj_Hough(Point p);

    float Get_Theta(){return theta_hough;};
    float Get_Phi(){return phi_hough;};
    
    WCVertex* get_vertex(){return vertex;};
    WCTrack* get_track(){return track;};

    bool IsShower(MergeSpaceCellSelection& mcells);
    bool IsContained(WCShower *shower);
    bool Contain(WCTrack *track);
    bool Overlap(WCShower *shower);
    


  protected:
    WCVertex *vertex;
    WCTrack *track;
    MergeSpaceCellSelection all_mcells;
    MergeSpaceCellSelection& exclude_mcells;
    MergeSpaceCellMap& mcells_map;

    float theta_hough;
    float phi_hough;

    int ncell_out_range;
    int ncell_angle_out_range;
    int ncell_total;
    int nmcell_out_range;
    int nmcell_total;
  };
  
  typedef std::vector<WCShower*> WCShowerSelection;
  
}

#endif
