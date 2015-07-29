#ifndef ClusterTrack_h
#define ClusterTrack_h

#include "WireCellData/MergeSpaceCell.h"
#include "TH2F.h"
#include <vector>
#include <map>


namespace WireCell {
  
  class ClusterTrack {
  public:
    ClusterTrack(MergeSpaceCell *cell);
    ~ClusterTrack();//{};

    void AddMSCell(MergeSpaceCell *cell);
    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};

    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

    void SC_Hough(Point& p, float dis = -1);
    void SC_Hough(Point& p1, Point&p, float dis = -1);
    
    void SC_IterativeHough(Point &p, float dis = 3 * units::cm);

    bool CrossAll(Point &p, float theta, float phi);

    float Get_Theta();
    float Get_Phi();


  protected:
    MergeSpaceCellSelection all_mcells;
    

    TH2F *hough;
    /* std::vector<double> sc_theta; */
    /* std::vector<double> sc_phi; */
    /* std::vector<double> sc_q; */

  };
  
  typedef std::vector<ClusterTrack*> ClusterTrackSelection;
  typedef std::map<ClusterTrack*,MergeSpaceCellSelection> ClusterMSpaceCellMap;
  typedef std::map<MergeSpaceCell*,ClusterTrackSelection> MSpaceCellClusterMap;
}

#endif
