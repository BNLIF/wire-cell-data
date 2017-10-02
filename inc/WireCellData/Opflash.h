#ifndef Opflash_h
#define Opflash_h

#include "WireCellData/COphit.h"
#include "TH1F.h"

namespace WireCell{
  class Opflash{
  public:
    Opflash(COphitSelection &ophits);
    Opflash(TH1F **hist, double start_time, int start_bin, int end_bin, float bin_width=6*15.625/1000.);
	    
    ~Opflash();

    double get_time(){return time;};
    double get_total_PE(){return total_PE;};
    double get_PE(int ch){return PE[ch];};
    double get_PE_err(int ch){return PE_err[ch];};
    bool get_fired(int ch);
    int get_num_fired(){return fired_channels.size();};

    int get_type(){return type;}
    double get_low_time(){return low_time;};
    double get_high_time(){return high_time;};
    
  protected:
    
    int type;
    double low_time;
    double high_time;
    
    double time;
    double total_PE;
    std::vector<int> fired_channels;
    double PE[32];
    double PE_err[32];
  };
  
  typedef std::vector<Opflash*> OpflashSelection;
}

#endif
