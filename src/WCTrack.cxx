#include "WireCellData/WCTrack.h"

using namespace WireCell;

WCTrack::WCTrack(MergeClusterTrack& mct)
  : mct(mct)
{
}

WCTrack::~WCTrack(){
}


int WCTrack::TrackType(MergeSpaceCell& cell){
  int result = 0;
  return result;
}
