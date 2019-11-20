#include "WCPData/Opflash.h"
#include <algorithm>
#include <iostream>

using namespace WCP;

WCP::Opflash::Opflash(int type, int flash_id, double low_time, double high_time, double time, double total_PE, std::vector<int>& fired_channels, double* temp_PE, double* temp_PE_err, std::vector<double>& l1_fired_time, std::vector<double>& l1_fired_pe)
  : type(type)
  , flash_id(flash_id)
  , low_time(low_time)
  , high_time(high_time)
  , time(time)
  , total_PE(total_PE)
  , fired_channels(fired_channels)
  , l1_fired_time(l1_fired_time)
  , l1_fired_pe(l1_fired_pe)
{
  for (int i=0;i!=32;i++){
    PE[i] = *(temp_PE+i);
    PE_err[i] = *(temp_PE_err+i);
  }
}

WCP::Opflash::Opflash(COphitSelection &ophits)
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

  // std::cout << total_PE << std::endl;
  if (total_PE !=0){
    time /= total_PE;
  }else{
    time = ophits.at(0)->get_time();
  }

  low_time = time - 3 * 15.625/1000.;
  high_time = time + 37 * 15.625/1000.;
  
}

WCP::Opflash::Opflash(TH1F **hist, double start_time, int start_bin, int end_bin, float bin_width)
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
  int max_bin = start_bin;

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
  double pe_noise_scaled = 1.875*(end_bin-std::max(start_bin,0))/78.;
  for (int i=0;i!=32;i++){
    PE[i] -= pe_noise_scaled;
    if (PE[i]<0) PE[i] = 0;
    total_PE += PE[i];

    if (PE[i]>=2){
      fired_channels.push_back(i);
    }

    PE_err[i] = sqrt(pow(PE_err[i],2) // standard error
		     + (PE[i] + pe_noise_scaled*2) //statistical term
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

WCP::Opflash::~Opflash(){
  
}

void WCP::Opflash::swap_channels(){
  //harded coded for now ... 
  for (auto it = fired_channels.begin(); it!=fired_channels.end(); it++){
    if (*it == 26) *it = 27;
    if (*it == 27) *it = 28;
    if (*it == 28) *it = 29;
    if (*it == 29) *it = 30;
    if (*it == 30) *it = 31;
    if (*it == 31) *it = 26;
  }
  
  double temp;
  temp = PE[31];
  PE[31] = PE[30];
  PE[30] = PE[29];
  PE[29] = PE[28];
  PE[28] = PE[27];
  PE[27] = PE[26];
  PE[26] = temp;

  temp = PE_err[31];
  PE_err[31] = PE_err[30];
  PE_err[30] = PE_err[29];
  PE_err[29] = PE_err[28];
  PE_err[28] = PE_err[27];
  PE_err[27] = PE_err[26];
  PE_err[26] = temp;
}

void WCP::Opflash::Add_l1info(TH1F *hist_tot_pe, TH1F *hist_mult, double start_time, int start_bin, int end_bin, float bin_width){
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


bool WCP::Opflash::get_fired(int ch){
  if (std::find(fired_channels.begin(),fired_channels.end(),ch)==fired_channels.end()){
    return false;
  }else{
    return true;
  }
}
