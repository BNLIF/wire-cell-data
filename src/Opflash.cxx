#include "WireCellData/Opflash.h"

using namespace WireCell;

WireCell::Opflash::Opflash(COphitSelection &ophits){
  for (int i=0;i!=32;i++){
    PE[i] = 0;
    PE_err[i] = 4.6;
  }
  time = 0;
  total_PE = 0;
  for (size_t i=0; i!= ophits.size(); i++){
    fired_channels.push_back(ophits.at(i)->get_ch_no());
    PE[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE();
    PE_err[ophits.at(i)->get_ch_no()] = ophits.at(i)->get_PE_err();
    time += ophits.at(i)->get_PE() * ophits.at(i)->get_time();
    total_PE += ophits.at(i)->get_PE() ;
  }
  time /= total_PE;
}

WireCell::Opflash::~Opflash(){
  
}

bool WireCell::Opflash::get_fired(int ch){
  if (find(fired_channels.begin(),fired_channels.end(),ch)==fired_channels.end()){
    return false;
  }else{
    return true;
  }
}
