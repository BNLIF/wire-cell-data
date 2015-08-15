#ifndef WCVertex_h
#define WCVertex_h

#include "WireCellData/Point.h"
#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/WCTrack.h"

//#include "Minuit2/Minuit2Minimizer.h"
//#include "Math/Functor.h"
#include "Minuit2/FCNBase.h"



#include <vector>
#include <map>

namespace WireCell {
  class WCVertex {
    
  public:
    Point Center();
    WCVertex(MergeSpaceCell& msc);
    ~WCVertex();
    void Add(WCTrack* track);
    int get_ntracks(){return tracks.size();};
    WCTrackSelection& get_tracks(){return tracks;};
    MergeSpaceCell* get_msc(){return msc;};

    int IsInside(WCVertex *vertex);
    bool AddVertex(WCVertex *vertex, int flag = 1);
    bool CheckContain(MergeSpaceCell *cell);
    void OrganizeTracks();
    WCTrackSelection BreakTracks();
    void ProcessTracks(WCTrackSelection& break_tracks);
    void OrganizeEnds(MergeSpaceCellSelection& cells, int flag = 1);
    // static double dis2(const double *xx);
    int FindVertex();
    
    double get_ky(int i){return tracks_ky.at(i);}
    double get_kz(int i){return tracks_kz.at(i);}
    

  protected:
    Point center;
    MergeSpaceCell *msc;
    WCTrackSelection tracks;
    std::vector<double> tracks_ky;
    std::vector<double> tracks_kz;
  };
  
  typedef std::vector<WCVertex*> WCVertexSelection;
  typedef std::map<MergeSpaceCell*, WCVertex*> MSC_WCV_Map;
  //typedef std::map<WCTrack*, WCVertexSelection> WCT_WCVs_Map;


  class MyFCN : public ROOT::Minuit2::FCNBase { 
    
  public: 
    double Up() const { return 1.; }
    
    MyFCN(WCVertex* vertex) 
      : vertex(vertex)
      {}
    
    ~MyFCN(){}
        
    double operator() (const std::vector<double> & xx) const;
    double get_chi2(const std::vector<double> & xx) const;

  private:
    WCVertex *vertex;
  };

}

#endif
