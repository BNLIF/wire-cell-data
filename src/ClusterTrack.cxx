#include "WireCellData/ClusterTrack.h"

using namespace WireCell;

ClusterTrack::ClusterTrack(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
}

ClusterTrack::~ClusterTrack(){
}


void ClusterTrack::AddMSCell(MergeSpaceCell *cell){
  all_mcells.push_back(cell);
}


