#include "WireCellData/Slim3DDeadCluster.h"

using namespace std;
using namespace WireCell;

Slim3DDeadCluster::Slim3DDeadCluster(SlimMergeGeomCell &cell, int time_slice){
  
}

Slim3DDeadCluster::~Slim3DDeadCluster(){
  cluster.clear();
  gcluster.clear();
}

int Slim3DDeadCluster::AddCell(SlimMergeGeomCell &cell, int time_slice, int offset){
  
}


void Slim3DDeadCluster::MergeCluster(Slim3DDeadCluster &cluster1){

}
