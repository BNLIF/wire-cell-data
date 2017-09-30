#ifndef Opflash_h
#define Opflash_h

#include "WireCellData/COphit.h"

namespace WireCell{
  class Opflash{
  public:
    Opflash(COphitSelection &ophits);
    ~Opflash();

    double get_time(){return time;};
    double get_total_PE(){return total_PE;};
    double get_PE(int ch){return PE[ch];};
    double get_PE_err(int ch){return PE_err[ch];};
    bool get_fired(int ch);
    int get_num_fired(){return fired_channels.size();};
    
  protected:
    double time;
    double total_PE;
    std::vector<int> fired_channels;
    double PE[32];
    double PE_err[32];
  };
  
  typedef std::vector<Opflash*> OpflashSelection;
}

#endif
