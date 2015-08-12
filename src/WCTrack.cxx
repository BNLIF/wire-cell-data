#include "WireCellData/WCTrack.h"

using namespace WireCell;

WCTrack::WCTrack(MergeClusterTrack& mct)
  : mct(mct)
{
  MergeSpaceCellSelection& mcells = mct.Get_allmcells();
  end_scells.push_back(mcells.front());
  
  MergeSpaceCell* last_cell = mct.Get_LastMSCell();

  for (int i = mct.Get_allmcells().size()-1;i>=0; i--){
    MergeSpaceCell *cell = mct.Get_allmcells().at(i);
    if (cell->Get_Center().x!=last_cell->Get_Center().x){
      end_scells.push_back(mct.Get_allmcells().at(i+1));
      break;
    }
  }

}

WCTrack::~WCTrack(){
}


int WCTrack::TrackType(MergeSpaceCell& cell){
  int result = 0;
  return result;
}
