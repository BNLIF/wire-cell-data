#ifndef WCTrack_h
#define WCTrack_h

#include "WireCellData/MergeClusterTrack.h"
//#include "WireCellData/WCVertex.h"
#include <vector>
#include <map>

namespace WireCell {
  class WCTrack{
  public:
    WCTrack(MergeClusterTrack& mct);
    ~WCTrack();
    MergeClusterTrack& get_mct(){return mct;};
    int TrackType(MergeSpaceCell& cell);
    MergeSpaceCellSelection& get_end_scells(){return end_scells;};
    MergeSpaceCellSelection& get_all_cells(){return all_cells;};

    MergeSpaceCell* replace_end_scells(MergeSpaceCell* cell, MergeSpaceCellSelection* cells=0);
    void ModifyCells();
    void ReplaceEndCell(MergeSpaceCell *cell1, MergeSpaceCell *cell2);
    bool Grow(MergeSpaceCell *cell, int flag = 1); 

  protected:
    MergeClusterTrack& mct;
    MergeSpaceCellSelection end_scells;
    MergeSpaceCellSelection all_cells;
    

  };

  typedef std::vector<WCTrack*> WCTrackSelection;
  typedef std::map<MergeClusterTrack*,WCTrack*> MCT_WCT_Map;
  //typedef std::map<WCTrack*, WCVertexSelection> WCT_WCV_Map;
}

#endif
