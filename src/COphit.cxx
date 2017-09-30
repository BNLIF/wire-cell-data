#include "WireCellData/COphit.h"

using namespace WireCell;

WireCell::COphit::COphit(int ch_no, TH1S *hist, double time, double gain)
  : channel_no(ch_no)
  , time(time)
  , gain(gain)
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

  
}

WireCell::COphit::~COphit(){
}
