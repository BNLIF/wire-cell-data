#include "WireCellData/Projected2DCluster.h"

using namespace WireCell;

Projected2DCluster::Projected2DCluster(WirePlaneType_t plane_no, GeomCell *parent_cell)
  : plane_no(plane_no)
  , parent_cell(parent_cell)
{
}

Projected2DCluster::~Projected2DCluster(){
}


void Projected2DCluster::AddCell(SlimMergeGeomCell *mcell){
  
}

