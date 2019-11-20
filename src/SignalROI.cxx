#include "WCPData/SignalROI.h"

#include <algorithm>

using namespace WCP;

SignalROI::SignalROI(WirePlaneType_t plane, int chid, int start_bin, int end_bin, TH1F *h1)
  : plane(plane)
  , chid(chid)
  , start_bin(start_bin)
  , end_bin(end_bin)
{
  float start_content = h1->GetBinContent(start_bin+1);
  float end_content = h1->GetBinContent(end_bin+1);
  contents.resize(end_bin-start_bin+1);

  
  for (int i=start_bin; i<= end_bin; i++){
    float content = h1->GetBinContent(i+1) - ((end_content - start_content)*(i-start_bin)/(end_bin-start_bin) + start_content);
    contents.at(i-start_bin) = content;
  }
}

SignalROI::SignalROI(SignalROI *roi){
  plane = roi->get_plane();
  chid = roi->get_chid();
  start_bin = roi->get_start_bin();
  end_bin = roi->get_end_bin();
  for (int i=start_bin; i<=end_bin;i++){
    contents.push_back(roi->get_contents().at(i-start_bin));
  }
}

double SignalROI::get_average_heights(){
  double sum1 = 0;
  double sum2 = 0;
  for (size_t i=0;i!=contents.size();i++){
    sum1 += contents.at(i);
    sum2 ++;
  }
  if (sum2!=0){
    return sum1/sum2;
  }else{
    return 0;
  }
}


bool SignalROI::overlap(SignalROI* roi1, float th, float th1){
  int min_start_bin = start_bin;
  if (start_bin < roi1->get_start_bin())
    min_start_bin = roi1->get_start_bin();
  int min_end_bin = end_bin;
  if (end_bin > roi1->get_end_bin())
    min_end_bin = roi1->get_end_bin();
  if (min_end_bin > min_start_bin){
    std::vector<float>& contents1 = roi1->get_contents();

    for (int i=min_start_bin; i<= min_end_bin; i++){
      if (contents.at(i-start_bin) > th && 
	  contents1.at(i-roi1->get_start_bin())>th1){
	  return true;
      }
    }
    
    return false;
  }else{
    return false;
  }
}

bool SignalROI::overlap(SignalROI* roi){

  int min_start_bin = start_bin;
  if (start_bin < roi->get_start_bin())
    min_start_bin = roi->get_start_bin();
  int min_end_bin = end_bin;
  if (end_bin > roi->get_end_bin())
    min_end_bin = roi->get_end_bin();
  if (min_end_bin > min_start_bin){
    return true;
  }else{
    return false;
  }
  // if (start_bin < roi->get_end_bin() && end_bin >= roi->get_end_bin())
  //   return true;
  // if (start_bin <= roi->get_start_bin() && end_bin > roi->get_start_bin())
  //   return true;
  // if (start_bin >= roi->get_start_bin() && end_bin <= roi->get_end_bin())
  //   return true;
  //return false;
}



SignalROI::~SignalROI(){
}



std::vector<std::pair<int,int>> SignalROI::get_above_threshold(float th){
  std::vector<std::pair<int,int>> bins;
  for (int i=0;i<contents.size();i++){
    if (contents.at(i) > th){
      int start = i;
      int end = i;
      for (int j=i+1;j<contents.size();j++){
	if (contents.at(j) > th){
	  end = j;
	}else{
	  break;
	}
      }
      bins.push_back(std::make_pair(start,end));
      i = end;
    }
  }

  return bins;
}
