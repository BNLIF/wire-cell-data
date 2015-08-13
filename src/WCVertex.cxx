#include "WireCellData/WCVertex.h"

using namespace WireCell;

WCVertex::WCVertex(MergeSpaceCell& msc)
  : msc(msc)
{
  center = msc.Get_Center();
}

WCVertex::~WCVertex(){
}


Point WCVertex::Center(){
  return center;
}

void WCVertex::Add(WCTrack* track){
  tracks.push_back(track);
}


int WCVertex::IsInside(WCVertex *vertex){
  int result = 1;
  WCTrackSelection& temp_tracks = vertex->get_tracks();
  for (int i=0;i!=tracks.size();i++){
    WCTrack *track = tracks.at(i);
    auto it = find(temp_tracks.begin(),temp_tracks.end(),track);
    if (it == temp_tracks.end()){
      result = -1;
      break;
    }
  }
  
  if (result == 1 && tracks.size() == temp_tracks.size()){
    result = 0;
  }
  

  return result;
}
