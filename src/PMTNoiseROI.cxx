#include "WireCellData/PMTNoiseROI.h"

#include <algorithm>

using namespace WireCell;


PMTNoiseROI::PMTNoiseROI(int start_bin, int end_bin)
  : start_bin(start_bin)
  , end_bin(end_bin)
{
  
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
  

  return true;
}
