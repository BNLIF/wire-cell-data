#ifndef WIRECELL_LMBDT_H
#define WIRECELL_LMBDT_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TH1D.h"

namespace WireCell{
  class LMBDT{
  public:
    LMBDT(TH1D* signalEff, TH1D* backgroundEff, double bdtScore);
    ~LMBDT();

    void set_BDT_score(double value);
    
    double get_BDT_score_eff_signal(double value);
    double get_BDT_score_eff_background(double value);
    double get_BDT_score_max_significance(); // S/sqrt(S+B)

    void callReader(double *pred_PE, double *flash_PE, double *max_PE, double *ks_dis, double *chi2, double *ndf, double *cluster_length,
		    int *runNo, int *subRunNo, int *eventNo, int *flash_id, int *tpc_cluster_id, int *flag_anode, int *flag_boundary,
		    int *trigger_type, int *scan_type, int *event_type);

    bool isSignal();
    
  private:
    TH1D* signalEff;
    TH1D* backgroundEff;
    double bdtScore;
  };
}

#endif
  
