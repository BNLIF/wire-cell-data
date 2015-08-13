#ifndef WCVertex_h
#define WCVertex_h

#include "WireCellData/Point.h"
#include "WireCellData/MergeSpaceCell.h"
#include "WireCellData/WCTrack.h"

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
    MergeSpaceCell* get_msc(){return &msc;};

    int IsInside(WCVertex *vertex);
    bool AddVertex(WCVertex *vertex);
    
    bool CheckContain(MergeSpaceCell *cell);

    void OrganizeTracks();

  protected:
    Point center;
    MergeSpaceCell& msc;
    WCTrackSelection tracks;
  };
  
  typedef std::vector<WCVertex*> WCVertexSelection;
  typedef std::map<MergeSpaceCell*, WCVertex*> MSC_WCV_Map;
  //typedef std::map<WCTrack*, WCVertexSelection> WCT_WCVs_Map;
}

#endif
