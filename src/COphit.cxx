#include "WireCellData/COphit.h"

using namespace WireCell;

WireCell::COphit::COphit(int ch_no, TH1S *hist, double time, double gain, double gain_err)
  : channel_no(ch_no)
  , time(time)
  , gain(gain)
  , gain_err(gain_err)
{
  // calculate baseline
  baseline = hist->GetBinContent(1);
  
  // calculate peak and integral 
  peak = 0;
  integral = 0;
  for (int i=0; i!=40; i++){
    double content = hist->GetBinContent(i+1) - baseline;
    if (content > peak){
      peak  = content;
    }
    integral += content;
  }

  // calculate PE and its error ... 
  good_baseline = false;
  if (integral > 20e3 || fabs(baseline-2050)<50) good_baseline = true;

  // special treatment of FEM channel 28 ... 
  if (channel_no == 28) integral = cal_integral(peak);
  if (!good_baseline) integral = cal_integral(peak);
  
  if (good_baseline){
    PE = integral / gain * 2; // taking into account factor of 2 for 0.6 us window ... 
    PE_err = sqrt(pow(PE * gain_err/gain,2) + pow(PE*0.02,2)); // standard error below threshold would be 4.6 PE ... 8/sqrt(3)

    if (PE < 12){ // too small ... 
      // some issues ... 
      PE_err = sqrt(pow(PE_err,2) + pow(4.6,2));
      good_baseline = false;
    }
    
  }else{
    PE = integral/gain * 2;
    PE_err = sqrt(pow(PE,2)+4.6*4.6); // given 100% uncertainties ... 
  }
  
}

WireCell::COphit::~COphit(){
}

double WireCell::COphit::cal_integral(double peak){
  double content;
  if (peak <=4000){
    content = 1.06241e+01 * peak + 2.01214e-04 * pow(peak,2);
  }else{
    content = 4.5715824e4 + 8.62296e+00 * (peak - 4000) + 6.76898e-04 * pow(peak-4000,2);
  }
  return content;
}
