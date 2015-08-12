#ifndef WCTrack_h
#define WCTrack_h

#include "WireCellData/MergeClusterTrack.h"
#include "WireCellData/WCVertex.h"
#include <vector>
#include <map>

namespace WireCell {
  class WCTrack{
  public:
    WCTrack(MergeClusterTrack& mct);
    ~WCTrack();
    MergeClusterTrack& get_mct(){return mct;};
    
    int TrackType(MergeSpaceCell& cell);
    
  protected:
    MergeClusterTrack& mct;
    
  };

  typedef std::vector<WCTrack*> WCTrackSelection;
  typedef std::map<MergeClusterTrack*,WCTrack*> MCT_WCT_Map;
  //typedef std::map<WCTrack*, WCVertexSelection> WCT_WCV_Map;
}

#endif
