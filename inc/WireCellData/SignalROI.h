#ifndef SignalROI_h
#define SignalROI_h

#include "WireCellData/GeomWire.h"

#include "TH1F.h"

#include <iostream>
#include <vector>
#include <list>
#include <map>


namespace WireCell{
  class SignalROI{
   public:
    SignalROI(WirePlaneType_t plane, int chid, int start_bin, int end_bin, TH1F *h1);
    SignalROI(SignalROI *roi);
    ~SignalROI();
    int get_start_bin(){return start_bin;}
    int get_end_bin(){return end_bin;}
    int get_chid(){return chid;}
    WirePlaneType_t get_plane(){return plane;}
    std::vector<float>& get_contents(){return contents;}
    std::vector<std::pair<int,int>> get_above_threshold(float th);
  
    bool overlap(SignalROI *roi);
    
  private:
    int start_bin;
    int end_bin;
    int chid;
    WirePlaneType_t plane;
    std::vector<float> contents;
  };
  
  typedef std::list<SignalROI*>SignalROIList;
  typedef std::vector<SignalROI*> SignalROISelection; 
  typedef std::vector<SignalROISelection> SignalROIChSelection;
  typedef std::vector<SignalROIList> SignalROIChList;
  typedef std::map<SignalROI*, SignalROISelection> SignalROIMap;
  
}

#endif
