#ifndef WCShower_h
#define WCShower_h

#include "WCPData/WCVertex.h"
#include "WCPData/MergeSpaceCell.h"

namespace WCP{
  
  class WCShower{
  public:
    WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cells,MergeSpaceCellMap& mcells_map);
    ~WCShower();
	    
    MergeSpaceCellSelection& get_all_cells(){return all_mcells;};

    void Iterate(MergeSpaceCell *curr_cell, WCP::MergeSpaceCellSelection &curr_cells);
    
    void SC_Hough(Point p);
    float SC_proj_Hough(Point p);

    float Get_Theta(){return theta_hough;};
    float Get_Phi(){return phi_hough;};
    
    WCVertex* get_vertex(){return vertex;};
    WCTrack* get_track(){return track;};

    bool IsShower(MergeSpaceCellSelection& mcells);
    bool IsShower_CheckVertex();
    bool IsContained(WCShower *shower);
    bool Contain(WCTrack *track);
    bool Overlap(WCShower *shower);
    
    float distance_to_center();

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
