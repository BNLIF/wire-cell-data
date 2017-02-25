#include "WireCellData/PMTNoiseROI.h"

#include <algorithm>

using namespace WireCell;


PMTNoiseROI::PMTNoiseROI(int start_bin, int end_bin, int peak)
  : start_bin(start_bin)
  , end_bin(end_bin)
{
  peaks.push_back(peak);
}


PMTNoiseROI::~PMTNoiseROI(){
  
}

void PMTNoiseROI::insert_peak(int peak){
  if (find(peaks.begin(),peaks.end(),peak)==peaks.end())
    peaks.push_back(peak);
}


void PMTNoiseROI::insert_uwires(int wire_no){
  if (find(induction_uwires.begin(),induction_uwires.end(),wire_no)==induction_uwires.end()){
    induction_uwires.push_back(wire_no);
  }
}

void PMTNoiseROI::insert_vwires(int wire_no){
  if (find(induction_vwires.begin(),induction_vwires.end(),wire_no)==induction_vwires.end()){
    induction_vwires.push_back(wire_no);
  }
}


bool PMTNoiseROI::merge_ROI(PMTNoiseROI& ROI){
  // decide how to merge two ROIs? ... 

  // one peak is contained in the other's range then merge
  if (peaks.at(0)>=ROI.get_start_bin() && peaks.at(0)<=ROI.get_end_bin()
      || start_bin <= ROI.get_peaks().at(0) && end_bin >= ROI.get_peaks().at(0)){
    if (ROI.get_start_bin() < start_bin)
      start_bin = ROI.get_start_bin();
    if (ROI.get_end_bin() > end_bin)
      end_bin = ROI.get_end_bin();
    
    for (int i=0;i!=ROI.get_peaks().size();i++){
      if (find(peaks.begin(),peaks.end(),ROI.get_peaks().at(i)) == peaks.end())
	peaks.push_back(ROI.get_peaks().at(i));
    }

    
    return true;
  }
 
  return false;
}
