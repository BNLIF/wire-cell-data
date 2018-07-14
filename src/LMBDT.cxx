#include "WireCellData/LMBDT.h"

using namespace WireCell;
using namespace TMVA;

LMBDT::LMBDT(float pred_PE, float flash_PE, float max_PE, float ks_dis,
	     float chi2, float ndf, float cluster_length,
	     int event_type, int flag_anode, int flag_boundary)
  : bdtScore(-1)
{
  reader = new TMVA::Reader( "!Color:!Silent" );
  reader->AddVariable("mag := log(pred_PE/flash_PE)",   &mag);
  reader->AddVariable("shape := ks_dis",                &shape);
  reader->AddVariable("chi := log(chi2/ndf)",           &chi);
  reader->AddVariable("max := log(max_PE/flash_PE)",    &max);
  reader->AddVariable("cluster := log(cluster_length)", &cluster);
  reader->AddVariable("pred := pred_PE",                &pred);

  reader->AddSpectator("event_type",    &event_type);
  reader->AddSpectator("flag_anode",    &flag_anode);
  reader->AddSpectator("flag_boundary", &flag_boundary);

  TString method = "BDT";
  TString prefix = "input_data_files/lmBDT";
  TString methodName = method + TString(" method");
  TString weightFile = prefix + TString("_") + method + TString(".weights.xml");
  reader->BookMVA( methodName, weightFile );

  mag = log(pred_PE/flash_PE);
  shape = ks_dis;
  chi = log(chi2/ndf);
  max = log(max_PE/flash_PE);
  cluster = log(cluster_length);
  pred = pred_PE;
}

LMBDT::~LMBDT(){
  delete reader;
}

void LMBDT::set_BDT_score(double value){
  bdtScore = value;
}

double LMBDT::get_BDT_score_eff_signal(double value, TH1D* signalEff){

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

double LMBDT::get_BDT_score_eff_background(double value, TH1D* backgroundEff){

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

double LMBDT::get_BDT_score_max_significance(TH1D* signalEff, TH1D* backgroundEff){

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

bool LMBDT::isSignal(){

  if(reader->EvaluateMVA( "BDT method" )>bdtScore) return true;
  
  return false;
}

bool LMBDT::isSignal(double value){

  if(reader->EvaluateMVA( "BDT method" )>value) return true;

  return false;
}

