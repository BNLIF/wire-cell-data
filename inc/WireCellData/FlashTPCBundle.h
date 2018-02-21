#ifndef WIRECELL_FLASHTPCBUNDLE_H
#define WIRECELL_FLASHTPCBUNDLE_H

#include "WireCellData/PR3DCluster.h"
#include "WireCellData/Opflash.h"

namespace WireCell{
  class FlashTPCBundle{
  public:
    FlashTPCBundle(Opflash* flash, PR3DCluster *main_cluster, int flash_index_id, int cluster_index_id);
    ~FlashTPCBundle();

    void set_flag_close_to_PMT(bool value){flag_close_to_PMT = value;};
    void set_flag_at_x_boundary(bool value){flag_at_x_boundary = value;};
    
    bool get_flag_close_to_PMT(){return flag_close_to_PMT;};
    bool get_flag_at_x_boundary(){return flag_at_x_boundary;};

    std::vector<double>& get_pred_pmt_light(){return pred_pmt_light;};
    Opflash* get_flash(){return flash;};
    PR3DCluster* get_main_cluster(){return main_cluster;};
    PR3DClusterSelection& get_other_clusters(){return other_clusters;};
    PR3DClusterSelection& get_more_clusters(){return more_clusters;};
    
  private:
    Opflash *flash;
    PR3DCluster *main_cluster;

    int cluster_index_id;
    int flash_index_id;

    bool flag_close_to_PMT;
    bool flag_at_x_boundary;
    
    std::vector<double> pred_pmt_light; // prediction

    PR3DClusterSelection other_clusters; // save every other one 
    PR3DClusterSelection more_clusters;  // save ones satisfying the cut    
    
    
  };

  typedef std::vector<FlashTPCBundle*> FlashTPCBundleSelection;
  typedef std::map<Opflash*, FlashTPCBundleSelection> Flash_bundles_map;
  typedef std::map<PR3DCluster*, FlashTPCBundleSelection> Cluster_bundles_map;
  
}


#endif
