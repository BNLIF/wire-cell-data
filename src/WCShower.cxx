#include "WireCellData/WCShower.h"

using namespace WireCell;

WCShower::WCShower(WCVertex *vertex, WCTrack *track, MergeSpaceCellSelection& exclude_cell,MergeSpaceCellMap& mcells_map)
  : vertex(vertex)
  , track(track)
{
  
}

WCShower::~WCShower(){
}
