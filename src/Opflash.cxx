#include "WireCellData/Opflash.h"
#include <algorithm>
#include <iostream>

using namespace WireCell;

WireCell::Opflash::Opflash(COphitSelection &ophits)
  : type(1)
  , flash_id (-1)
{
  for (int i=0;i!=32;i++){
    PE[i] = 0;
    PE_err[i] = 6.4; // 11/sqrt(3.) fudge now ... 
  }
  time = 0;
  total_PE = 0;
  for (size_t i=0; i!= ophits.size(); i++){
    fired_channels.push_back(ophits.at(i)->get_ch_no());
    PE[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE() - 0.15*2; // 250 kHz at 0.6 us ... 
    PE_err[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE_err();
    time += ophits.at(i)->get_PE() * ophits.at(i)->get_time();
    total_PE += ophits.at(i)->get_PE() ;
  }
  time /= total_PE;

  low_time = time - 3 * 15.625/1000.;
  high_time = time + 37 * 15.625/1000.;
  
}

WireCell::Opflash::Opflash(TH1F **hist, double start_time, int start_bin, int end_bin, float bin_width)
  : type(2)
  , flash_id (-1)
{
  low_time = start_time + (start_bin+0.5)*bin_width;
  high_time = start_time + (end_bin-0.5)*bin_width;
  
  for (int i=0;i!=32;i++){
    PE[i] = 0;
    PE_err[i] = 0.2; // 0.2 PE as base ... 
  }
  
  double max = 0;
  int max_bin = 0;

  for (int i=start_bin; i!=end_bin;i++){
    double peak = 0;
    double mult = 0;
    for (int j=0;j!=32;j++){
      double content = hist[j]->GetBinContent(i+1);
      if (content < 0.2) content = 0;
      peak += content;
      PE[j] += content;
      if (content>1.5) {
	//	std::cout << i << " " << j << " " << content << std::endl;
	mult++;
      }
    }

    if (peak > max && mult >=3){
      //std::cout << max << " " << peak << " " << max_bin << " " << i << " " << mult << std::endl;
      max = peak;
      max_bin = i;
    }
  }

  time = start_time + (max_bin+0.5) * bin_width;

  total_PE = 0;
  for (int i=0;i!=32;i++){
    PE[i] -= 1.875; // 7.5 us * random noise ... 
    if (PE[i]<0) PE[i] = 0;
    total_PE += PE[i];

    if (PE[i]>=2){
      fired_channels.push_back(i);
    }

    PE_err[i] = sqrt(pow(PE_err[i],2) // standard error
		     + (PE[i] + 1.875*2) //statistical term
		     + pow(PE[i]*0.02,2) // 2% base systematic uncertainties
		     );
    
    // if (PE[i] < 1000){
    //   PE_err[i] = sqrt(pow(PE_err[i],2) + pow(PE[i]*0.02,2)); // standard error below threshold would be 4.6 PE ... 8/sqrt(3)
    // }else if (PE[i] < 2000){
    //   PE_err[i] = sqrt(pow(PE_err[i],2) + pow(PE[i]*0.06,2)); // standard error below threshold would be 4.6 PE ... 8/sqrt(3)
    // }else{
    //   PE_err[i] = sqrt(pow(PE_err[i],2) + pow(PE[i]*0.18,2)); // standard error below threshold would be 4.6 PE ... 8/sqrt(3)
    // }
    
    
  }
  
}

WireCell::Opflash::~Opflash(){
  
}

void WireCell::Opflash::Add_l1info(TH1F *hist_tot_pe, TH1F *hist_mult, double start_time, int start_bin, int end_bin, float bin_width){
  std::vector<int> fired_bin;
  std::vector<double> fired_pe;
  for (int i=start_bin; i!=end_bin;i++){
    double pe = hist_tot_pe->GetBinContent(i+1);
    double mult = hist_mult->GetBinContent(i+1);
    if (pe >=10 && mult>=3){
      fired_bin.push_back(i);
      fired_pe.push_back(pe);
    }
  }

  for (size_t i=0; i!=fired_bin.size();i++){
    if (i==0){
      double time = start_time + (fired_bin.at(0)+0.5)*bin_width;
      l1_fired_time.push_back(time);
      l1_fired_pe.push_back(fired_pe.at(0));
    }else{
      if (fired_bin.at(i)-fired_bin.at(i-1)==1){
	l1_fired_pe.back() += fired_pe.at(i);
      }else{
	double time = start_time + (fired_bin.at(i)+0.5)*bin_width;
	l1_fired_time.push_back(time);
	l1_fired_pe.push_back(fired_pe.at(i));
      }
    }
  }

  // std::cout << l1_fired_time.size() << std::endl;
  // for (size_t i = 0; i!= l1_fired_time.size(); i++){
  //   std::cout << l1_fired_time.at(i) << " " << l1_fired_pe.at(i) << std::endl;
  // }
  
}


bool WireCell::Opflash::get_fired(int ch){
  if (std::find(fired_channels.begin(),fired_channels.end(),ch)==fired_channels.end()){
    return false;
  }else{
    return true;
  }
}
