#ifndef COphit_h
#define COphit_h

#include "TH1S.h"
#include <vector>


namespace WireCell{

  class COphit{
  public:
    COphit(int ch_no, TH1S *hist, double time, double gain, double gain_err);
    ~COphit();

    double get_time(){return time;};
    double get_baseline(){return baseline;};
    double get_peak(){return peak;};
    double get_integral(){return integral;};
    double get_gain(){return gain;}
    int get_ch_no(){return channel_no;};

    // derivated ... 
    double get_PE(){return PE;};
    double get_PE_err(){return PE_err;};


    
  protected:
    double cal_integral(double peak);

    
    int channel_no; // FEM channel number
    double gain; // to be set
    double gain_err;
    
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
