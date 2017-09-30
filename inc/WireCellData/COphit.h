#ifndef COphit_h
#define COphit_h

#include "TH1S.h"
#include <vector>


namespace WireCell{

  class COphit{
  public:
    COphit(int ch_no, TH1S *hist, double time, double gain);
    ~COphit();

    double get_time(){return time;};
    double get_baseline(){return baseline;};
    double get_peak(){return peak;};
    double get_integral(){return integral;};
    double get_gain(){return gain;}
    int get_ch_no(){return channel_no;};

   
    
  protected:
    int channel_no; // FEM channel number
    double gain; // to be set
    
    double time; // start time in us ... 
    double baseline; // baseline 
    double peak; // maximum PE
    double integral; // integral 

    double PE;
    double PE_err;
    
  };

  typedef std::vector<COphit*> COphitSelection;
  
}

#endif
