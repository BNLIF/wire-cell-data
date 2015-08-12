#include "WireCellData/WCVertex.h"

using namespace WireCell;

WCVertex::WCVertex(MergeSpaceCell& msc)
  : msc(msc)
{
}

WCVertex::~WCVertex(){
}


Point WCVertex::Center(){
  return center;
}
