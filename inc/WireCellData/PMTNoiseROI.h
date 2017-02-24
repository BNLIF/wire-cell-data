#ifndef PMTNoiseROI_h
#define PMTNoiseROI_h

#include <vector>
#include <list>

namespace WireCell{
  class PMTNoiseROI{
  public:
    PMTNoiseROI(int start_bin, int end_bin);
    ~PMTNoiseROI();

    int get_start_bin(){return start_bin;};
    int get_end_bin(){return end_bin;};
    
    void insert_peak(int peak);
    
    void insert_uwires(int wire_no);
    void insert_vwires(int wire_no);

    std::vector<int>& get_peaks(){return peaks;}
    std::list<int>& get_uwires(){return induction_uwires;}
    std::list<int>& get_vwires(){return induction_vwires;}

    bool merge_ROI(PMTNoiseROI& ROI);

  private: 
    std::vector<int> peaks;
    int start_bin;
    int end_bin;
    
    std::list<int> induction_uwires;
    std::list<int> induction_vwires;
  };
}

#endif
