#include "WCPData/TPCParams.h"
#include "TFile.h"

#include <fstream>
#include <iostream>

using namespace WCP;

  // set defaults
WCP::TPCParams::TPCParams()
  : m_pitch_u(3)
  , m_pitch_v(3)
  , m_pitch_w(3)
  , m_ts_width(3.2)
  , m_angle_u(1.0472)
  , m_angle_v(-1.0472)
  , m_angle_w(0)
  , first_u_dis(0)
  , first_v_dis(0)
  , first_w_dis(0)
  , nrebin(4) 
  , time_offset(4)
  , flag_corr(false)
  , gu(0)
  , gv(0)
  , gw(0)
  , g_proton(0)
  , g_muon(0)
  , g_pion(0)
  , g_kaon(0)
  , g_electron(0)
{
  //init_corr_files();
};

WCP::TPCParams::~TPCParams(){
  if (gu!=0) delete gu;
  if (gv!=0) delete gv;
  if (gw!=0) delete gw;

  if (g_proton!=0) delete g_proton;
  if (g_electron!=0) delete g_electron;
  if (g_muon!=0) delete g_muon;
  if (g_pion!=0) delete g_pion;
  if (g_kaon!=0) delete g_kaon;
}

void WCP::TPCParams::init_PID_dq_dx(TString filename){
  if (g_proton!=0) delete g_proton;
  if (g_electron!=0) delete g_electron;
  if (g_muon!=0) delete g_muon;
  if (g_pion!=0) delete g_pion;
  if (g_kaon!=0) delete g_kaon;
  
  TFile *file = new TFile("./input_data_files/stopping_ave_dQ_dx.root");
  g_muon = (TGraph*)file->Get("muon");
  g_pion = (TGraph*)file->Get("pion");
  g_kaon = (TGraph*)file->Get("kaon");
  g_proton = (TGraph*)file->Get("proton");
  g_electron = (TGraph*)file->Get("electron");
  file->Close();
}

void WCP::TPCParams::init_corr_files(TString file_u , int ndata_u, TString file_v, int ndata_v, TString file_w, int ndata_w){
  std::cout << "Initialize dQ/dx Correction Files" << std::endl;
  if (gu!=0) delete gu; gu = new TGraph();
  if (gv!=0) delete gv; gv = new TGraph();
  if (gw!=0) delete gw; gw = new TGraph();
  double x,y;
  std::ifstream file1(file_u);
  for (int i=0;i!=ndata_u;i++){
    file1 >> x >> y;
    gu->SetPoint(i,x,y);
  }
  std::ifstream file2(file_v);
  for (int i=0;i!=ndata_v;i++){
    file2 >> x >> y;
    gv->SetPoint(i,x,y);
  }
  std::ifstream file3(file_w);
  for (int i=0;i!=ndata_w;i++){
    file3 >> x >> y;
    gw->SetPoint(i,x,y);
  }
  flag_corr = true;
}

double WCP::TPCParams::get_corr_factor(WCP::Point& p, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw){
  double u = offset_u + slope_yu * p.y + slope_zu * p.z;
  double v = offset_v + slope_yv * p.y + slope_zv * p.z;
  double w = offset_w + slope_yw * p.y + slope_zw * p.z;
  //  std::cout << u << " " << v << " " << w << std::endl;
  double corr = gu->Eval(u) * gv->Eval(v) * gw->Eval(w);
  // corr = 1;
  return corr;
  
}
