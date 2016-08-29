#ifndef SignalROI_h
#define SignalROI_h

#include "TH1F.h"

#include <iostream>
#include <vector>
#include <map>


namespace WireCell{
  class SignalROI{
   public:
    SignalROI(int start_bin, int end_bin, TH1F *h1);
    ~SignalROI();
    int get_start_bin(){return start_bin;}
    int get_end_bin(){return end_bin;}
    std::vector<float>& get_contents(){return contents;}
    std::vector<std::pair<int,int>> get_above_threshold(float th);
  
  private:
    int start_bin;
    int end_bin;
    std::vector<float> contents;
  };
  
  typedef std::vector<SignalROI*> SignalROISelection; 
  typedef std::map<SignalROI*, SignalROISelection> SignalROIMap;
}

#endif
