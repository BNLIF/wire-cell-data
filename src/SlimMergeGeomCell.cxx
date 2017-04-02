#include "WireCellData/SlimMergeGeomCell.h"

using namespace WireCell;

WireCell::SlimMergeGeomCell::SlimMergeGeomCell(){
  
}

WireCell::SlimMergeGeomCell::~SlimMergeGeomCell(){
  uwires.clear();
  vwires.clear();
  wwires.clear();
}

void WireCell::SlimMergeGeomCell::AddWire(GeomWire *wire, WirePlaneType_t plane){
  if (plane == WirePlaneType_t(0)){
    if (find(uwires.begin(),uwires.end(),wire)==uwires.end())
      uwires.push_back(wire);
  }else if (plane == WirePlaneType_t(1)){
    if (find(vwires.begin(),vwires.end(),wire)==vwires.end())
      vwires.push_back(wire);
  }else if (plane == WirePlaneType_t(2)){
    if (find(wwires.begin(),wwires.end(),wire)==wwires.end())
      wwires.push_back(wire);
  }
}


void WireCell::SlimMergeGeomCell::OrderWires(){
  WireCell::sort_by_ident(uwires);
  WireCell::sort_by_ident(vwires);
  WireCell::sort_by_ident(wwires);
}
