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
