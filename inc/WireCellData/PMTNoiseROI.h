#ifndef PMTNoiseROI_h
#define PMTNoiseROI_h

#include <vector>

namespace WireCell{
  class PMTNoiseROI{
  public:
    PMTNoiseROI(int start_bin, int end_bin, int peak);
    ~PMTNoiseROI();

    int get_start_bin(){return start_bin;};
    int get_end_bin(){return end_bin;};
    
    void insert_peak(int peak);
    
    void insert_uwires(int wire_no);
    void insert_vwires(int wire_no);

    std::vector<int>& get_peaks(){return peaks;}
    
    std::vector<int>& get_uwires(){return induction_uwires;}
    std::vector<int>& get_vwires(){return induction_vwires;}
    
    std::vector<int>& get_sorted_uwires(){return sorted_ind_uwires;}    
    std::vector<int>& get_sorted_vwires(){return sorted_ind_vwires;}    
    
    void sort_wires(int nwire=1);
    bool merge_ROI(PMTNoiseROI& ROI);

  private: 
    std::vector<int> peaks;
    int start_bin;
    int end_bin;
    
    std::vector<int> induction_uwires;
    std::vector<int> induction_vwires;
    
    std::vector<int> sorted_ind_uwires;
    std::vector<int> sorted_ind_vwires;
  };
}

#endif
