#include "WCPData/PhotonLibrary.h"

#include "TGraph.h"
#include "TChain.h"

#include <iostream>

WCP::Photon_Library::Photon_Library(double eventTime, int run_no, bool flag_data, bool flag_add_light_yield_err, bool flag_timestamp){


  //  std::cout << "ZXin_1: " << eventTime << " " << flag_timestamp << std::endl;
  
  rel_light_yield_err = 0;
  scaling_light_mag = 0.01 * 1.5;
  Double_t yield_run_no[37]={5590, 5934, 6207, 6427, 6617, 6854, 7059, 7305, 7648, 8199, 8518, 8871, 9209, 9468, 9652, 10478, 10701, 10924, 11197, 11605, 11816, 12021, 12344, 12505, 13521, 13725, 14034, 14256, 14527, 14773, 15013, 15426, 15922, 16218, 16643, 16977, 17417};
  Double_t yield_time[37]={1458959696,1460824531,1462698133,1464589035,1466452319,
			   1468322746,1470354984,1472048832,1473667444,1476001396,
			   1477674527,1479544558,1481410591,1483293418,1485007076,
			   1488911371,1490773998,1492647736,1494220283,1496560957,
			   1498242134,1499863075,1502616255,1503261543,1508234775,
			   1509464364,1511351145,1513205887,1515074419,1516950264,
			   1518816572,1520690576,1522560834,1524425522,1526295309,
			   1528168585,1530046184};
  Double_t yield_ratio[37]={1, 1, 1, 0.99, 1.01, 0.97, 0.97, 0.97, 0.97, 0.96, 0.97, 0.9, 0.83, 0.82, 0.8, 0.77, 0.77, 0.77, 0.76, 0.71, 0.73, 0.71, 0.7, 0.68, 0.64, 0.65, 0.65, 0.64, 0.64, 0.63, 0.63, 0.64, 0.64, 0.64, 0.64, 0.62, 0.62};
  Double_t yield_ratio_err[37]={0, 0, 0, 0.01, 0, 0.01, 0.01, 0, 0.01, 0.02, 0.01, 0.03, 0.05, 0.06, 0.06, 0.07, 0.05, 0.08, 0.08, 0.09, 0.09, 0.1, 0.16, 0.1, 0.1, 0.1, 0.11, 0.1, 0.11, 0.11, 0.11, 0.1, 0.1, 0.11, 0.11, 0.1, 0.1};
  Int_t start_run_no = 5590;
  Int_t end_run_no = 17417;
  
  TGraph gratio;//(37, yield_run_no, yield_ratio);
  TGraph gratio_err;//(37, yield_run_no, yield_ratio_err);
  for (int i=0;i!=37;i++){
    if (flag_timestamp){
      gratio.SetPoint(i,yield_time[i],yield_ratio[i]);
      gratio_err.SetPoint(i,yield_time[i],yield_ratio_err[i]);
    }else{
      gratio.SetPoint(i,yield_run_no[i],yield_ratio[i]);
      gratio_err.SetPoint(i,yield_run_no[i],yield_ratio_err[i]);
    }
  }
  
//  std::map<int,int> map_lib_pmt,map_pmt_lib;
  map_lib_pmt[1]=2; map_pmt_lib[2]=1; 
  map_lib_pmt[0]=4; map_pmt_lib[4]=0; 
  map_lib_pmt[3]=0; map_pmt_lib[0]=3; 
  map_lib_pmt[2]=5; map_pmt_lib[5]=2; 
  map_lib_pmt[5]=1; map_pmt_lib[1]=5; 
  map_lib_pmt[4]=6; map_pmt_lib[6]=4; 
  map_lib_pmt[6]=3; map_pmt_lib[3]=6; 
  
  map_lib_pmt[9]=7; map_pmt_lib[7]=9; 
  map_lib_pmt[7]=9; map_pmt_lib[9]=7; 
  map_lib_pmt[8]=11; map_pmt_lib[11]=8; 
  map_lib_pmt[11]=8; map_pmt_lib[8]=11;  
  map_lib_pmt[10]=12; map_pmt_lib[12]=10; 
  map_lib_pmt[12]=10; map_pmt_lib[10]=12; 

  map_lib_pmt[14]=13; map_pmt_lib[13]=14;  
  map_lib_pmt[13]=15; map_pmt_lib[15]=13; 
  map_lib_pmt[15]=17; map_pmt_lib[17]=15; 
  map_lib_pmt[17]=14; map_pmt_lib[14]=17; 
  map_lib_pmt[16]=18; map_pmt_lib[18]=16; 
  map_lib_pmt[18]=16; map_pmt_lib[16]=18; 

  map_lib_pmt[21]=19; map_pmt_lib[19]=21; 
  map_lib_pmt[22]=20; map_pmt_lib[20]=22; 
  map_lib_pmt[19]=21; map_pmt_lib[21]=19; 
  map_lib_pmt[20]=23; map_pmt_lib[23]=20; 
  map_lib_pmt[23]=24; map_pmt_lib[24]=23; 
  map_lib_pmt[24]=22; map_pmt_lib[22]=24; 

  map_lib_pmt[26]=25; map_pmt_lib[25]=26; 
  map_lib_pmt[27]=30; map_pmt_lib[30]=27; 
  map_lib_pmt[28]=31; map_pmt_lib[31]=28; 
  map_lib_pmt[31]=29; map_pmt_lib[29]=31;
  // original map
  map_lib_pmt[25]=28; map_pmt_lib[28]=25; 
  map_lib_pmt[30]=27; map_pmt_lib[27]=30; 
  map_lib_pmt[29]=26; map_pmt_lib[26]=29;

  // fixed map ... (if not swap in the flash reconstruction ...)
  // map_lib_pmt[25]=27; map_pmt_lib[27]=25; 
  // map_lib_pmt[30]=26; map_pmt_lib[26]=30; 
  // map_lib_pmt[29]=28; map_pmt_lib[28]=29;
  
  TChain *T = new TChain("/pmtresponse/PhotonLibraryData","/pmtresponse/PhotonLibraryData");
  T->AddFile("./uboone_photon_library.root");
  //std::cout << T->GetEntries();
  if (T->GetEntries()<=0) {
      std::cout << "Error: 2dToy::tpc_light_match failed to read uboone_photon_library.root"<<std::endl;
      exit(1);
  }
  Int_t Voxel;
  Int_t OpChannel;
  Float_t Visibility;
  T->SetBranchAddress("Voxel",&Voxel);
  T->SetBranchAddress("OpChannel",&OpChannel);
  T->SetBranchAddress("Visibility",&Visibility);

//  std::vector<std::list<std::pair<int,float>>> photon_library;
  library.resize(400*75*75);
  for (int i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    library.at(Voxel).push_back(std::make_pair(OpChannel,Visibility));
  }
  
  // hack
  //run_no = 17394;
  
  // for data, scale down the light ...
  if (flag_data){
    double scaling_additional = 1;
    double scaling_additional_err = 0;

    if (flag_timestamp){
      if (eventTime < yield_time[0]){
	scaling_additional = yield_ratio[0];
	scaling_additional_err = yield_ratio_err[0];
      }else if (eventTime > yield_time[36]){
	scaling_additional = yield_ratio[36];
	scaling_additional_err = yield_ratio_err[36];
      }else{
	scaling_additional = gratio.Eval(eventTime);
	scaling_additional_err = gratio_err.Eval(eventTime);
      }
    }else{
      if (run_no < start_run_no){
	scaling_additional = yield_ratio[0];
	scaling_additional_err = yield_ratio_err[0];
      }else if (run_no > end_run_no){
	scaling_additional = yield_ratio[36];
	scaling_additional_err = yield_ratio_err[36];
      }else{
	scaling_additional = gratio.Eval(run_no);
	scaling_additional_err = gratio_err.Eval(run_no);
      }
    }

    //std::cout << "kaka: " << scaling_additional << " " << scaling_additional_err << " " << flag_timestamp << std::endl;
    
    scaling_light_mag *= scaling_additional;
    if (flag_add_light_yield_err)
      rel_light_yield_err = scaling_additional_err/scaling_additional;
  }
}
