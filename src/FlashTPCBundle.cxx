#include "WireCellData/FlashTPCBundle.h"

using namespace WireCell;

FlashTPCBundle::FlashTPCBundle(Opflash* flash, PR3DCluster *main_cluster, int flash_index_id, int cluster_index_id)
  : flash(flash)
  , main_cluster(main_cluster)
  , flash_index_id(flash_index_id)
  , cluster_index_id(cluster_index_id)
  , flag_close_to_PMT(false)
  , flag_at_x_boundary(false)
{
  pred_pmt_light.resize(32,0);
}

FlashTPCBundle::~FlashTPCBundle(){
  
}

