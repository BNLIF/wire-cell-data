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
    void set_msc(MergeSpaceCell * cell){msc=cell;};
    void set_center(Point &p){
      center.x = p.x;
      center.y = p.y;
      center.z = p.z;
    }

    int IsInside(WCVertex *vertex);
    bool AddVertex(WCVertex *vertex, int flag = 1);
    bool CheckContain(MergeSpaceCell *cell);
    void OrganizeTracks();
    WCTrackSelection BreakTracks();
    WCTrackSelection BreakTracksAngle(WCTrackSelection& finished_tracks);

    void ProcessTracks(WCTrackSelection& break_tracks);
    void OrganizeEnds(MergeSpaceCellSelection& cells, int flag = 1);
    // static double dis2(const double *xx);
    bool FindVertex(int flag = 0);

   
    double get_ky(WCTrack *track);
    double get_kz(WCTrack *track);

    double get_ky(int i){return tracks_ky.at(i);}
    double get_kz(int i){return tracks_kz.at(i);}
    int get_fit_type(){return fit_type;};
    bool get_fit_success(){return fit_success;};
    bool get_fit_chi2(){return fit_chi2;};
    void reset_center();
    
  protected:
    int fit_type;

    bool fit_success;
    double fit_chi2;

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
    
  MyFCN(WCVertex* vertex, int fit_flag) 
      : vertex(vertex)
      , fit_flag(fit_flag)
      {}
    
    ~MyFCN(){}
        
    double operator() (const std::vector<double> & xx) const;
    double get_chi2(const std::vector<double> & xx) const;

  private:
    WCVertex *vertex;
    int fit_flag;
  };

}

#endif
