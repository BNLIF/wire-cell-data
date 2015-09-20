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

    bool AddMSCell(MergeSpaceCell *cell);
    MergeSpaceCellSelection& Get_allmcells(){return all_mcells;};

    MergeSpaceCell* Get_FirstMSCell(){return all_mcells.front();};
    MergeSpaceCell* Get_LastMSCell(){return all_mcells.back();};

    void SC_Hough(Point& p, float dis = -1, int flag = 1);
    void SC_Hough(Point& p1, Point&p, float dis = -1, int flag = 1);
    Point SC_IterativeHough(Point &p, float dis = 3 * units::cm, int flag = 1);

    int CrossNum(Point &p, float theta, float phi);
    int CrossNum(MergeSpaceCell *mcell1, float theta, float phi);

    float Get_Theta();
    float Get_Phi();

    float Get_CTheta(Point &p);
    float Get_CPhi(Point &p);

  protected:
    MergeSpaceCellSelection all_mcells;
   
    float theta_hough;
    float phi_hough;

    //    TH2F *hough;
    /* std::vector<double> sc_theta; */
    /* std::vector<double> sc_phi; */
    /* std::vector<double> sc_q; */

  };
  
  typedef std::vector<ClusterTrack*> ClusterTrackSelection;
  typedef std::map<ClusterTrack*,MergeSpaceCellSelection> ClusterMSpaceCellMap;
  typedef std::map<MergeSpaceCell*,ClusterTrackSelection> MSpaceCellClusterMap;
}

#endif
