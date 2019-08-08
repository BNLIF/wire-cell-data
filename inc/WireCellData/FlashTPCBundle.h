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
    void set_flash(Opflash* flash1){flash=flash1;};
    double get_total_pred_light();
    
    PR3DCluster* get_main_cluster(){return main_cluster;};
    void set_main_cluster(PR3DCluster* cluster){main_cluster = cluster;};
    PR3DClusterSelection& get_other_clusters(){return other_clusters;};
    PR3DClusterSelection& get_more_clusters(){return more_clusters;};

    bool examine_bundle(Double_t *cos_pe_low, Double_t *cos_pe_mid);
    bool examine_bundle(FlashTPCBundle* bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid);
    bool examine_bundle_rank(FlashTPCBundle* bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid);
    void add_bundle(FlashTPCBundle* bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid);
    bool examine_beam_bundle();

    //  bool check_tgm( WireCell2dToy::ToyFiducial *fid, double offset_x );
    
    
    double get_chi2(){return chi2;};
    int get_ndf(){return ndf;};
    double get_ks_dis(){return ks_dis;};
    void set_consistent_flag(bool value){flag_high_consistent = value;};
    bool get_consistent_flag(){return flag_high_consistent;};

    void set_spec_end_flag(bool value){flag_spec_end = value;};
    bool get_spec_end_flag(){return flag_spec_end;};
    bool get_potential_bad_match_flag(){return flag_potential_bad_match;};

    double get_strength(){return strength;};
    void set_strength(double value){ strength = value;};

    void examine_merge_clusters(double dis_cut = 3.6*units::cm);
    
  private:
    Opflash *flash;
    PR3DCluster *main_cluster;

    double strength;
    
    int cluster_index_id;
    int flash_index_id;

    bool flag_close_to_PMT;
    bool flag_at_x_boundary;
    bool flag_spec_end;

    bool flag_potential_bad_match;
    
    std::vector<double> pred_pmt_light; // prediction

    double ks_dis;
    double chi2;
    int ndf;
    bool flag_high_consistent;

    // add some varialbes for LM events ...

    // place holder for Brooke ...
    
    PR3DClusterSelection other_clusters; // save every other one 
    PR3DClusterSelection more_clusters;  // save ones satisfying the cut    
  };

  typedef std::vector<FlashTPCBundle*> FlashTPCBundleSelection;

  struct FlashTPCBundleComparePtr {
	bool operator() (FlashTPCBundle* a, FlashTPCBundle* b) {
	  if (a && b) {
	    if (a->get_flash() && b->get_flash()){
	      if (a->get_flash()->get_flash_id() < b->get_flash()->get_flash_id()){
		return true;
	      }else if (a->get_flash()->get_flash_id() == b->get_flash()->get_flash_id()){
		if (a->get_main_cluster()->get_cluster_id() < b->get_main_cluster()->get_cluster_id()){
		  return true;
		}else{
		  return false;
		}
	      }else{
		return false;
	      }
	    }else if (a->get_flash()==0 && b->get_flash()){
	      return false;
	    }else if (a->get_flash() && b->get_flash()==0){
	      return true;
	    }else if (a->get_flash()==0 && b->get_flash()==0){
	      if (a->get_main_cluster()->get_cluster_id() < b->get_main_cluster()->get_cluster_id()){
		return true;
	      }else{
		return false;
	      }
	    }else{
	      return false;
	    }
	  }
	  return false;
	}
    };

  struct OpflashComparePtr {
    bool operator() (Opflash* a, Opflash* b) {
      if (a && b) {
	if (a->get_flash_id() < b->get_flash_id()){
	  return true;
	}else{
	  return false;
	}
      }
      return false;
    }
  };

   struct PR3DClusterComparePtr {
    bool operator() (PR3DCluster* a, PR3DCluster* b) {
      if (a && b) {
	if (a->get_cluster_id() < b->get_cluster_id()){
	  return true;
	}else{
	  return false;
	}
      }
      return false;
    }
  };

  
  
  
  typedef std::set<FlashTPCBundle*, FlashTPCBundleComparePtr> FlashTPCBundleSet;
  typedef std::map<Opflash*, FlashTPCBundleSelection, OpflashComparePtr > Flash_bundles_map;
  typedef std::map<PR3DCluster*, FlashTPCBundleSelection, PR3DClusterComparePtr > Cluster_bundles_map;
  
}


#endif
