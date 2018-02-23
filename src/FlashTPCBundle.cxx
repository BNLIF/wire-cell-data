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

void FlashTPCBundle::examine_bundle(Double_t *cos_pe_low, Double_t *cos_pe_mid){
  std::cout << flash->get_type() << " " << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << std::endl;
  for (int i=0;i!=32;i++){
    std::cout << flash->get_PE(i) << " ";
  }
  std::cout << std::endl;
  for (int i=0;i!=32;i++){
    std::cout << flash->get_PE_err(i) << " ";
  }
  std::cout << std::endl;
  for (int i=0;i!=32;i++){
    std::cout << pred_pmt_light.at(i) << " ";
  }
  std::cout << std::endl;
}
