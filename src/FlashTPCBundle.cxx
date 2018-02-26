#include "WireCellData/FlashTPCBundle.h"

using namespace WireCell;

FlashTPCBundle::FlashTPCBundle(Opflash* flash, PR3DCluster *main_cluster, int flash_index_id, int cluster_index_id)
  : flash(flash)
  , main_cluster(main_cluster)
  , flash_index_id(flash_index_id)
  , cluster_index_id(cluster_index_id)
  , flag_close_to_PMT(false)
  , flag_at_x_boundary(false)
  , ks_dis(1)
  , chi2(0)
  , ndf(0)
  , flag_high_consistent(false)
{
  pred_pmt_light.resize(32,0);
}

FlashTPCBundle::~FlashTPCBundle(){
  
}

bool FlashTPCBundle::examine_bundle(FlashTPCBundle *bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid){
  TH1F *h1 = new TH1F("h1","h1",32,0,32);
  TH1F *h2 = new TH1F("h2","h2",32,0,32);

  double pe[32],pe_err[32];
  double pred_pe[32];
  
  
  for (int i=0;i!=32;i++){
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_pmt_light.at(i) + bundle->get_pred_pmt_light().at(i);
  }
  
  for (int j=0;j!=32;j++){
    h1->SetBinContent(j+1,pe[j]);
    if ((pred_pe[j] < cos_pe_low[j] ||
	 (pred_pe[j] < cos_pe_mid[j]*1.1 && pe[j] ==0 )) && flash->get_type()==1){
      pred_pe[j] = 0;
    }
    h2->SetBinContent(j+1,pred_pe[j]);
  }


  double temp_chi2 = 0;
  double temp_ndf = 0;
  double temp_ks_dis = 1;
  
  double max_chi2 = 0;
  int max_bin = -1;
  
  if (h2->GetSum()!=0){
    temp_ks_dis = h1->KolmogorovTest(h2,"M");
  }
  
  for (int j=0;j!=32;j++){
    double cur_chi2 = 0;
    if (flag_close_to_PMT){
      if (pe[j]-pred_pe[j]>350&&pe[j]>pred_pe[j]*1.3){ // if the measurement is much larger than the prediction
	cur_chi2 = pow(pred_pe[j]-pe[j],2)/(pow(pe_err[j],2)+pow(pe[j]*0.5,2));
      }else{
	cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2); 
      }
    }else{
      cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2); 
    }
    temp_chi2 += cur_chi2;

    if (cur_chi2 > max_chi2){
      max_chi2 = cur_chi2;
      max_bin = j;
    }      
    
     if (pe[j]==0&&pred_pe[j]==0){
     }else{
       temp_ndf++;
     }
  }
  if (pe[max_bin] == 0 && pred_pe[max_bin]>0) // allow one PMT to be inefficient in measurement ... 
    temp_chi2 -= max_chi2-1;
  
  
  delete h1;
  delete h2;

  // if (main_cluster->get_cluster_id()==2)
  //   std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << ks_dis << " " << temp_ks_dis << " " << chi2 << " " << temp_chi2 << " " << ndf << std::endl;
  
  if ((temp_ks_dis < ks_dis + 0.06 &&
      temp_ks_dis < ks_dis * 1.2 && 
      temp_chi2 < chi2 + ndf * 5 &&
       temp_chi2 < chi2 * 1.21) ||
      (temp_ks_dis < ks_dis &&
       temp_chi2 < chi2 + ndf * 10 &&
       temp_chi2 < chi2 * 1.45)
      ){
    return true;
  }else{
    return false;
  }
  
}

bool FlashTPCBundle::examine_bundle(Double_t *cos_pe_low, Double_t *cos_pe_mid){
  TH1F *h1 = new TH1F("h1","h1",32,0,32);
  TH1F *h2 = new TH1F("h2","h2",32,0,32);

  double pe[32],pe_err[32];
  double pred_pe[32];

  for (int i=0;i!=32;i++){
    pe[i] = flash->get_PE(i);
    pe_err[i] = flash->get_PE_err(i);
    pred_pe[i] = pred_pmt_light.at(i);
  }

  for (int j=0;j!=32;j++){
    h1->SetBinContent(j+1,pe[j]);
    if ((pred_pe[j] < cos_pe_low[j] ||
	 (pred_pe[j] < cos_pe_mid[j]*1.1 && pe[j] ==0 )) && flash->get_type()==1){
      pred_pe[j] = 0;
    }
    h2->SetBinContent(j+1,pred_pe[j]);
  }

  if (h2->GetSum()!=0){
    ks_dis = h1->KolmogorovTest(h2,"M");
  }

  chi2 = 0;
  ndf = 0;
  double max_chi2 = 0;
  int max_bin = -1;
  
  for (int j=0;j!=32;j++){
    double cur_chi2 = 0;
    
    if (flag_close_to_PMT){
      if (pe[j]-pred_pe[j]>350&&pe[j]>pred_pe[j]*1.3){ // if the measurement is much larger than the prediction
	cur_chi2 = pow(pred_pe[j]-pe[j],2)/(pow(pe_err[j],2)+pow(pe[j]*0.5,2));
      }else{
	cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2); 
      }
    }else{
      cur_chi2 = pow(pred_pe[j]-pe[j],2)/pow(pe_err[j],2); 
    }
    chi2 += cur_chi2;

    if (cur_chi2 > max_chi2){
      max_chi2 = cur_chi2;
      max_bin = j;
    }      
    
     if (pe[j]==0&&pred_pe[j]==0){
     }else{
       ndf++;
     }
  }
  if (pe[max_bin] == 0 && pred_pe[max_bin]>0) // allow one PMT to be inefficient in measurement ... 
    chi2 -= max_chi2-1;
  

  if (ks_dis < 0.12 && ndf >=2 && chi2 < ndf * 25){
    flag_high_consistent = true;
  }else if (flag_at_x_boundary && ndf >=1 && chi2 < 9 * ndf && ks_dis < 0.12){
    flag_high_consistent = true;
  }else if (chi2 < 4 * ndf && ndf >=3 && ks_dis < 0.15){
    flag_high_consistent = true;
  }else if (ks_dis < 0.12 && ndf >=5 && chi2 < ndf * 55 && flag_close_to_PMT){
    flag_high_consistent = true;
  }
  
  // if (ks_dis<0.35)
  //   std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << flash->get_type() << " " << ks_dis << " " << chi2 << " " << ndf << std::endl;

  
  // for (int i=0;i!=32;i++){
  //   std::cout << flash->get_PE(i) << " ";
  // }
  // std::cout << std::endl;
  // for (int i=0;i!=32;i++){
  //   std::cout << flash->get_PE_err(i) << " ";
  // }
  // std::cout << std::endl;
  // for (int i=0;i!=32;i++){
  //   std::cout << pred_pmt_light.at(i) << " ";
  // }
  // std::cout << std::endl;

  delete h1;
  delete h2;

  if (flag_high_consistent){
    return true;
  }else{
    double ntot   = 0;
    double ntot1 = 0;
    double nfired = 0;
    double nfired1 = 0;
    for (int j=0;j!=32;j++){
      if (pred_pe[j] >0){
	ntot ++;
	if (pe[j] > 0.33*pred_pe[j])
	  nfired ++;
      }
      if (pred_pmt_light.at(j) > 1.0){
	ntot1 ++;
	if (pe[j] > 0.33*pred_pmt_light.at(j))
	  nfired1 ++;
      }
    }
    
    // if (fabs(main_cluster->get_cluster_id()-19)<=0 || main_cluster->get_cluster_id()==18){
    //std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << nfired << " " << ntot << " " << nfired1 << " " << ntot1 << " " << ks_dis << " " << chi2 << " " << ndf << " " << flag_at_x_boundary << " " << flag_close_to_PMT << std::endl;
    // }

    // exception for small clusters ... 
    if (nfired <=1 && (nfired1!=0 && nfired1 > 0.75*ntot1)) return true;
    
    if ( nfired==0 ) return false;
    if ( nfired < 0.5 * ntot && ntot - nfired >= 2) return false;
    if ( nfired==1 && ntot <=2 && nfired1< 0.2*ntot1) return false;
  }
  return true;
}
