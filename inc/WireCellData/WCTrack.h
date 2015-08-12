#ifndef WCTrack_h
#define WCTrack_h

#include "WireCellData/MergeClusterTrack.h"
#include <vector>

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
}

#endif
