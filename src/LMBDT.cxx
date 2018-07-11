#include "WireCellData/LMBDT.h"

//using namespace WireCell;
using namespace TMVA;

WireCell::LMBDT::LMBDT(TH1D* signalEff, TH1D* backgroundEff, double bdtScore)
  : signalEff(0)
  , backgroundEff(0)
  , bdtScore(-1)
{
}

WireCell::LMBDT::~LMBDT(){
}

void WireCell::LMBDT::set_BDT_score(double value){
  bdtScore = value;
}

double WireCell::LMBDT::get_BDT_score_eff_signal(double value){

  int lastBin = -1;
  int bins = signalEff->GetNbinsX();
  for(int i=1; i<=bins; i++){
    if(signalEff->GetBinContent(i)<value && i!=1){
      lastBin = i-1;
      break;
    }
  }
  
  bdtScore = signalEff->GetBinCenter(lastBin);
  return bdtScore;
}

double WireCell::LMBDT::get_BDT_score_eff_background(double value){

  int lastBin = -1;
  int bins = backgroundEff->GetNbinsX();
  for(int i=bins; i>=1; i--){
    if(backgroundEff->GetBinContent(i)>value && i!=bins){
      lastBin = i+1;
      break;
    }
  }
  
  bdtScore = backgroundEff->GetBinCenter(lastBin);
  return bdtScore;
}

double WireCell::LMBDT::get_BDT_score_max_significance(){

  int binsS = signalEff->GetNbinsX();
  int binsB = backgroundEff->GetNbinsX();
  if(binsS!=binsB){
    return bdtScore = -10;
  }

  double max = 0.;
  int lastBin = -1;
  for(int i=1; i<=binsS; i++){
    double temp = signalEff->GetBinContent(i) / sqrt(signalEff->GetBinContent(i)+backgroundEff->GetBinContent(i));
    if(temp > max && i!=1){
      max = temp;
      lastBin = i-1;
    }
    if(temp < max) break;
  }

  bdtScore = signalEff->GetBinCenter(lastBin);
  return bdtScore;
}

void WireCell::LMBDT::callReader(double *pred_PE, double *flash_PE, double *max_PE, double *ks_dis, double *chi2, double *ndf, double *cluster_length,
		       int *runNo, int *subRunNo, int *eventNo, int *flash_id, int *tpc_cluster_id, int *flag_anode, int *flag_boundary,
		       int *trigger_type, int *scan_type, int *event_type){

  //TMVA::Tools::Instance();
  
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  /*
  Float_t mag, shape, chi, max, cluster, pred;
  reader->AddVariable("mag := log(pred_PE/flash_PE)",   &mag);
  reader->AddVariable("shape := ks_dis",                &shape);
  reader->AddVariable("chi := log(chi2/ndf)",           &chi);
  reader->AddVariable("max := log(max_PE/flash_PE)",    &max);
  reader->AddVariable("cluster := log(cluster_length)", &cluster);
  reader->AddVariable("pred := pred_PE",                &pred);
  
  reader->AddSpectator("runNo", &runNo);
  reader->AddSpectator("subRunNo", &subRunNo);
  reader->AddSpectator("eventNo", &eventNo);
  reader->AddSpectator("flash_id", &flash_id);
  reader->AddSpectator("tpc_cluster_id", &tpc_cluster_id);
  reader->AddSpectator("scan_type", &scan_type);
  reader->AddSpectator("event_type", &event_type);
  reader->AddSpectator("trigger_type", &trigger_type);
  reader->AddSpectator("flag_anode", &flag_anode);
  reader->AddSpectator("flag_boundary", &flag_boundary);  
  */
}

bool WireCell::LMBDT::isSignal(){

  return false;
}

