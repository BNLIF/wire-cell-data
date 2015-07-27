#ifndef ClusterTrack_h
#define ClusterTrack_h

#include "WireCellData/MergeSpaceCell.h"
#include <vector>


namespace WireCell {
  
  class ClusterTrack {
  public:
    ClusterTrack(MergeSpaceCell *cell);
    ~ClusterTrack();//{};

    void AddMSCell(MergeSpaceCell *cell);
    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};

    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

  protected:
    MergeSpaceCellSelection all_mcells;

  };
  
  typedef std::vector<ClusterTrack*> ClusterTrackSelection;
}

#endif
