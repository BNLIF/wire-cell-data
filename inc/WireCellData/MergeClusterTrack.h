#ifndef MergeClusterTrack_h
#define MergeClusterTrack_h

#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/ClusterTrack.h"
#include "TH2F.h"
#include <vector>

namespace WireCell{
  class MergeClusterTrack{
  public:
    MergeClusterTrack(ClusterTrack *ctrack);
    ~MergeClusterTrack();
    
    void Add(ClusterTrack *ctrack, MergeSpaceCell *mcell1);

    
    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};
    void Update();

    void SC_Hough(Point& p, float dis = -1);
    void SC_Hough(Point& p1, Point&p, float dis = -1);
    
    Point SC_IterativeHough(Point &p, float dis = 3 * units::cm);

    ClusterTrack* GetClusterTrack(MergeSpaceCell* vertex);

    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

      float Get_Theta();
    float Get_Phi();


  protected:
    ClusterTrackSelection ctracks; // save merged clusters ... 
    MergeSpaceCellSelection all_mcells; // save all the merged cells
    MergeSpaceCellList all_mcells_list; // temporary one ... 
    
    TH2F *hough;
    
  };

  typedef std::vector<MergeClusterTrack*> MergeClusterTrackSelection;
  
}

#endif
