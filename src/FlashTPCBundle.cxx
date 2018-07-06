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
  , flag_spec_end(false)
  , flag_potential_bad_match(false)
  , strength(0)
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

  // if (main_cluster->get_cluster_id()==17)
  // std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << ks_dis << " " << temp_ks_dis << " " << chi2 << " " << temp_chi2 << " " << ndf << std::endl;
  
  if ((temp_ks_dis < ks_dis + 0.06 &&
       (temp_ks_dis < ks_dis * 1.2) && 
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

// bool FlashTPCBundle::check_tgm( WireCell2dToy::ToyFiducial *fid, double offset_x){
//   if (flash!=0){
//     // check the fiducial volume ...
//     std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
    
//     Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
//     Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
    
//     //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
//     bool flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
//     bool flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);
//     //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
    
//     // check the dead region ...
//     if (flag_inside_p1){
//       // define a local direction ...
//       TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
//       dir *= (-1);
//       flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
//     }
//     if (flag_inside_p2){
//       // define a  local direction ...
//       TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
//       dir *= (-1);
//       flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
//     }
    
//   }
  
//   return false;
// }

bool FlashTPCBundle::examine_bundle_rank(FlashTPCBundle *bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid){
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

  // if (main_cluster->get_cluster_id()==17)
  //  std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << bundle->get_main_cluster()->get_cluster_id() << " " << ks_dis << " " << temp_ks_dis << " " << chi2 << " " << temp_chi2 << " " << ndf << std::endl;
  
  if ((temp_ks_dis < ks_dis + 0.06 &&
       (temp_ks_dis < ks_dis * 1.2 || temp_ks_dis < 0.05 || temp_ks_dis < ks_dis + 0.03) && 
       temp_chi2 < chi2 + ndf * 5 &&
       temp_chi2 < chi2 * 1.21) ||
      (temp_ks_dis < ks_dis &&
       temp_chi2 < chi2 + ndf * 10 &&
       temp_chi2 < chi2 * 1.45) ||
      (temp_ks_dis * temp_chi2 < ks_dis * chi2)
      ){
    return true;
  }else{
    return false;
  }
  
}


void FlashTPCBundle::examine_merge_clusters(double dis_cut){
  // if (flash->get_type()==1) return; // only do for beam flash
  
  int main_cluster_id = main_cluster->get_cluster_id();

  PR3DClusterSelection merge_clusters;
  for (size_t i=0;i!=other_clusters.size();i++){
    PR3DCluster *temp_cluster = other_clusters.at(i);

    double dis_save = 1e9;
    
    {
      PR3DCluster *cluster1 = temp_cluster;
      PR3DCluster *cluster2 = main_cluster;
      SlimMergeGeomCell *prev_mcell1 = 0;
      SlimMergeGeomCell *prev_mcell2 = 0;
      SlimMergeGeomCell *mcell1 = 0; 
      Point p1;//
      SlimMergeGeomCell *mcell2=0;
      Point p2;
      
      mcell1 = *(cluster1->get_time_cells_set_map().begin()->second.begin());
      p1 = mcell1->center();
      
      while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
	prev_mcell1 = mcell1;
	prev_mcell2 = mcell2;
	
	// find the closest point and merged cell in cluster2
	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1);
	p2 = temp_results.second;
	mcell2 = temp_results.first;
	// find the closest point and merged cell in cluster1
	temp_results = cluster1->get_closest_point_mcell(p2);
	p1 = temp_results.second;
	mcell1 = temp_results.first;
      }
      double dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

      if (dis < dis_save){
	dis_save = dis;
      }

      prev_mcell1 = 0;
      prev_mcell2 = 0;
      
      mcell1 = *(cluster1->get_time_cells_set_map().rbegin()->second.begin());
      p1 = mcell1->center();

      while(mcell1!=prev_mcell1 || mcell2!=prev_mcell2){
	prev_mcell1 = mcell1;
	prev_mcell2 = mcell2;
	
	// find the closest point and merged cell in cluster2
	std::pair<SlimMergeGeomCell*,Point> temp_results = cluster2->get_closest_point_mcell(p1);
	p2 = temp_results.second;
	mcell2 = temp_results.first;
	// find the closest point and merged cell in cluster1
	temp_results = cluster1->get_closest_point_mcell(p2);
	p1 = temp_results.second;
	mcell1 = temp_results.first;
      }
      dis = sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));

      if (dis < dis_save){
	dis_save = dis;
      }
    }

    if (dis_save < dis_cut){
      merge_clusters.push_back(temp_cluster);
    }
    // if (main_cluster_id==13)
    //   std::cout << main_cluster_id << " " << temp_cluster->get_cluster_id() << " " << dis_save/units::cm << " " << merge_clusters.size() << std::endl;
  }

  if (merge_clusters.size()>0){
    merge_clusters.push_back(main_cluster);
    PR3DCluster *ncluster = new PR3DCluster(main_cluster_id);
    for (auto it1 = merge_clusters.begin(); it1!=merge_clusters.end(); it1++){
      PR3DCluster *ocluster = *(it1);
      SMGCSelection& mcells = ocluster->get_mcells();
      for (auto it2 = mcells.begin(); it2!=mcells.end(); it2++){
	SlimMergeGeomCell *mcell = (*it2);
	//std::cout << ocluster->get_cluster_id() << " " << mcell << std::endl;
	int time_slice = mcell->GetTimeSlice();
	ncluster->AddCell(mcell,time_slice);
      }
    }
    // delete old clusters
    for (auto it1 = merge_clusters.begin(); it1!=merge_clusters.end(); it1++){
      PR3DCluster *ocluster = *(it1);
      if (ocluster == main_cluster){
	delete ocluster;
      }else{
	auto it2 = find(other_clusters.begin(),other_clusters.end(),ocluster);
	if (it2!=other_clusters.end()){
	  other_clusters.erase(it2);
	}
	auto it3 = find(more_clusters.begin(), more_clusters.end(), ocluster);
	if (it3!=more_clusters.end()){
	  more_clusters.erase(it3);
	}
	delete ocluster;
      }
    }
    main_cluster = ncluster;
  }
  
  
}

void  FlashTPCBundle::add_bundle(FlashTPCBundle* bundle, Double_t *cos_pe_low, Double_t *cos_pe_mid){

  
  // flag_close_to_PMT = flag_close_to_PMT || bundle->get_flag_close_to_PMT();
  // flag_at_x_boundary = flag_at_x_boundary || bundle->get_flag_at_x_boundary();
  

  other_clusters.push_back(bundle->get_main_cluster());
  more_clusters.push_back(bundle->get_main_cluster());
  
  std::copy(bundle->get_other_clusters().begin(), bundle->get_other_clusters().end(), std::back_inserter(other_clusters));
  std::copy(bundle->get_more_clusters().begin(), bundle->get_more_clusters().end(), std::back_inserter(more_clusters));

  std::vector<double>& pes = bundle->get_pred_pmt_light();
  for (size_t i=0;i!=pred_pmt_light.size();i++){
    pred_pmt_light.at(i) += pes.at(i);
  }
  examine_bundle(cos_pe_low,cos_pe_mid);
}


bool FlashTPCBundle::examine_beam_bundle(){
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
    h2->SetBinContent(j+1,pred_pe[j]);
  }
  

  double temp_ks_dis = h1->KolmogorovTest(h2,"M");
  double temp_chi2 = 0;
  double temp_ndf = 0;
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
  h1->SetBinContent(max_bin+1,0);
  h2->SetBinContent(max_bin+1,0);
  double temp_ks_dis1 = h1->KolmogorovTest(h2,"M");

  //  std::cout << temp_ks_dis << " " << temp_ks_dis1 << " " << temp_chi2 << " " << temp_ndf << " " << max_chi2 << " " << max_bin << std::endl;

  if ( (temp_ks_dis < 0.1 ||temp_ks_dis1 < 0.05) &&
       (temp_chi2 < temp_ndf * 12 || temp_chi2 - max_chi2 < (temp_ndf-1)*6))
    return true;
  
  return false;
  
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
  
  if (ks_dis < 0.06 && ndf >=3 && chi2 < ndf * 36){
    flag_high_consistent = true;
  }else if (ks_dis<0.05 && ndf >=6 && chi2 < ndf * 45){
    flag_high_consistent = true;
  }else if (ks_dis < 0.12 && ndf >=3 && chi2 < ndf * 25){
    flag_high_consistent = true;
  }else if (flag_at_x_boundary && ndf >=2 && chi2 < 9 * ndf && ks_dis < 0.12){
    flag_high_consistent = true;
  }else if (flag_at_x_boundary && ndf >=1 && chi2 < 3 * ndf && ks_dis < 0.12){
    flag_high_consistent = true;
  }else if (chi2 < 4 * ndf && ndf >=3 && ks_dis < 0.15){
    flag_high_consistent = true;
  }else if (chi2 < 1.5 * ndf && ks_dis < 0.2 && ndf >= 3){
    flag_high_consistent = true;
  }else if (ks_dis < 0.12 && ndf >=5 && chi2 < ndf * 55 && flag_close_to_PMT){
    flag_high_consistent = true;
  }else if (ks_dis < 0.14 && ndf >=3 && chi2< ndf * 6){
    flag_high_consistent = true;
  }

 
  
  // if (ks_dis<0.35)
  // if (flash->get_flash_id()==14)
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
      if (pred_pe[j] >0.33){
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
    
    // if (fabs(main_cluster->get_cluster_id()-19)<=0){
    //   std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << nfired << " " << ntot << " " << nfired1 << " " << ntot1 << " " << ks_dis << " " << chi2 << " " << ndf << " " << flag_at_x_boundary << " " << flag_close_to_PMT << std::endl;
    // }

    // exception for small clusters ... 
    if (nfired <=1 && (nfired1!=0 && nfired1 > 0.75*ntot1)) return true;

    if (nfired <=2 && ks_dis>0.8 && chi2>60*ndf && ndf >=6) return false;
    

    if (flag_at_x_boundary && (!flag_close_to_PMT) ){
      if ( nfired==0 ) {
	flag_potential_bad_match = true;
      }
      if ( nfired==1 && ntot <=2 && nfired1< 0.2*ntot1) {
	flag_potential_bad_match = true;
      }
      if ( nfired < 0.5 * ntot && ntot - nfired >= 2) {
	flag_potential_bad_match = true;
      }
    }else{
      if ( nfired==0 ) {
	flag_potential_bad_match = true;
	return false;
      }
      if ( nfired==1 && ntot <=2 && nfired1< 0.2*ntot1) {
	flag_potential_bad_match = true;
	return false;
      }
      if ( nfired < 0.5 * ntot && ntot - nfired >= 2) {
	flag_potential_bad_match = true;
	return false;
      }
    }
  }
  return true;
}
