#ifndef MergeClusterTrack_h
#define MergeClusterTrack_h

#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/ClusterTrack.h"
#include <vector>

namespace WireCell{
  class MergeClusterTrack{
  public:
    MergeClusterTrack(ClusterTrack *ctrack);
    ~MergeClusterTrack();
    
    void Add(ClusterTrack *ctrack, MergeSpaceCell *mcell1);

    
    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};
    void Update();

    ClusterTrack* GetClusterTrack(MergeSpaceCell* vertex);

    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

  protected:
    ClusterTrackSelection ctracks; // save merged clusters ... 
    MergeSpaceCellSelection all_mcells; // save all the merged cells
    MergeSpaceCellList all_mcells_list; // temporary one ... 

  };

  typedef std::vector<MergeClusterTrack*> MergeClusterTrackSelection;
  
}

#endif
