#ifndef WIRECELL_LMBDT_H
#define WIRECELL_LMBDT_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TH1D.h"

namespace WireCell{
  class LMBDT{
  public:
    LMBDT(float pred_PE, float flash_PE, float max_PE, float ks_dis,
	  float chi2, float ndf, float cluster_length, 
	  int event_type, int flag_anode, int flag_boundary);
    ~LMBDT();

    void set_BDT_score(double value);
    
    double get_BDT_score_eff_signal(double value, TH1D* signalEff);
    double get_BDT_score_eff_background(double value, TH1D* backgroundEff);
    double get_BDT_score_max_significance(TH1D* signalEff, TH1D* backgroundEff); // S/sqrt(S+B)

    bool isSignal();
    bool isSignal(double value);
    
  private:
    
  protected:
    double bdtScore;
    
    TMVA::Reader *reader;
    
    float pred_PE;
    float flash_PE;
    float max_PE;
    float ks_dis;
    float chi2;
    float ndf;
    float cluster_length;

    int event_type;
    int flag_anode;
    int flag_boundary;

    float mag;
    float shape;
    float chi;
    float max;
    float cluster;
    float pred;
  };
}

#endif
  
