#ifndef WCTrack_h
#define WCTrack_h

#include "WireCellData/MergeClusterTrack.h"

namespace WireCell {
  class WCTrack{
  public:
    WCTrack(MergeClusterTrack& mct);
  protected:
    int track_type;
    MergeClusterTrack& mct;
  };
}

#endif
