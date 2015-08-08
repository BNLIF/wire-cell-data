#ifndef MergeClusterTrack_h
#define MergeClusterTrack_h

#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/ClusterTrack.h"
#include "TH2F.h"
#include <vector>
#include <map>
#include <set>

namespace WireCell{
  class MergeClusterTrack{
  public:
    MergeClusterTrack(ClusterTrack *ctrack);
    ~MergeClusterTrack();
    
    void Add(ClusterTrack *ctrack, MergeSpaceCell *mcell1);
    
    void Add(MergeSpaceCell *mcell, int flag);
    void Add(MergeClusterTrack *mctrack, MergeSpaceCell *mcell, int flag);
    
    void AddTrack(MergeClusterTrack *mctrack);
    void Organize();

    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};
    ClusterTrackSelection& Get_ctracks(){return ctracks;};
    
    void Update();

    void SC_Hough(Point& p, float dis = -1, int flag = 1);
    void SC_Hough(Point& p1, Point&p, float dis = -1, int flag = 1);
    Point SC_IterativeHough(Point &p, float dis = 3 * units::cm, int flag = 1);
    

    ClusterTrack* GetClusterTrack(MergeSpaceCell* vertex);
    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

    bool CheckCell(MergeSpaceCell* scell);

    float Get_Theta();
    float Get_Phi();

    int Get_TimeLength();
    MergeSpaceCellSelection& Get_MSCS(int time);
    int Get_Time(int time){return times.at(time);};

  protected:
    ClusterTrackSelection ctracks; // save merged clusters ... 
    MergeSpaceCellSelection all_mcells; // save all the merged cells
    MergeSpaceCellList all_mcells_list; // temporary one ... 

    std::vector<MergeSpaceCellSelection> times_mcells;
    std::vector<int> times;
    
    float theta_hough;
    float phi_hough;
    //TH2F *hough;
    
  };

  struct MergeClusterTrackCompare {
    bool operator() (MergeClusterTrack* a, MergeClusterTrack* b) const{
      if (a->Get_allmcells().size() == b->Get_allmcells().size()){
	return a > b;
      }
      return a->Get_allmcells().size() > b->Get_allmcells().size();
    }
  };

  typedef std::vector<MergeClusterTrack*> MergeClusterTrackSelection;
  typedef std::map<MergeSpaceCell*, MergeClusterTrackSelection> MSC_MCT_Map;
  typedef std::set<MergeClusterTrack*, MergeClusterTrackCompare> MergeClusterTrackSet;
}

#endif
