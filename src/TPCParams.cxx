#include "WCPData/TPCParams.h"
#include "WCPData/Units.h"
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
  , mass_proton(938.272*units::MeV)
  , mass_neutron(939.565*units::MeV)
  , mass_muon(105.658*units::MeV)
  , mass_pion(139.571*units::MeV)
  , mass_neutral_pion(134.977*units::MeV)
  , mass_kaon(493.677*units::MeV)
  , mass_electron(0.511*units::MeV)
  , g_proton_r2ke(0)
  , g_muon_r2ke(0)
  , g_pion_r2ke(0)
  , g_kaon_r2ke(0)
  , g_electron_r2ke(0)
  , electron_lifetime(1000)
  , h3_Dx(0)
  , h3_Dy(0)
  , h3_Dz(0)
  , h3_Ex(0)
  , h3_Ey(0)
  , h3_Ez(0)    
  , flag_PosEfield_corr(0)
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
  
  if (g_proton_r2ke!=0) delete g_proton_r2ke;
  if (g_electron_r2ke!=0) delete g_electron_r2ke;
  if (g_muon_r2ke!=0) delete g_muon_r2ke;
  if (g_pion_r2ke!=0) delete g_pion_r2ke;
  if (g_kaon_r2ke!=0) delete g_kaon_r2ke;
}

double WCP::TPCParams::get_attenuation_ratio(double drift_time){
  double ratio = 1;
  if (drift_time < 0) drift_time = 0;
  if (electron_lifetime >= 1000){
    return ratio;
  }else{
    //    ratio = exp(-drift_time/electron_lifetime + drift_time/1000.);
    ratio = exp(-drift_time/electron_lifetime );
    return ratio;
  }
}

void WCP::TPCParams::init_PID_dq_dx(TString filename, TString filename1){
  if (g_proton!=0) delete g_proton;
  if (g_electron!=0) delete g_electron;
  if (g_muon!=0) delete g_muon;
  if (g_pion!=0) delete g_pion;
  if (g_kaon!=0) delete g_kaon;

  if (g_proton_r2ke!=0) delete g_proton_r2ke;
  if (g_electron_r2ke!=0) delete g_electron_r2ke;
  if (g_muon_r2ke!=0) delete g_muon_r2ke;
  if (g_pion_r2ke!=0) delete g_pion_r2ke;
  if (g_kaon_r2ke!=0) delete g_kaon_r2ke;
  
  TFile *file = new TFile(filename);
  g_muon = (TGraph*)file->Get("muon");
  g_pion = (TGraph*)file->Get("pion");
  g_kaon = (TGraph*)file->Get("kaon");
  g_proton = (TGraph*)file->Get("proton");
  g_electron = (TGraph*)file->Get("electron");
  file->Close();

  TFile *file1 = new TFile(filename1);
  g_muon_r2ke = (TGraph*)file1->Get("muon");
  g_pion_r2ke = (TGraph*)file1->Get("pion");
  g_kaon_r2ke = (TGraph*)file1->Get("kaon");
  g_proton_r2ke = (TGraph*)file1->Get("proton");
  g_electron_r2ke = (TGraph*)file1->Get("electron");
  file1->Close();
}

void WCP::TPCParams::init_Pos_Efield_SCE_correction(TString filename)
{
  TFile *file = new TFile(filename, "read");
  TH3F *hDx = (TH3F*)file->Get("hDx");
  TH3F *hDy = (TH3F*)file->Get("hDy");
  TH3F *hDz = (TH3F*)file->Get("hDz");
  TH3F *hEx = (TH3F*)file->Get("hEx");
  TH3F *hEy = (TH3F*)file->Get("hEy");
  TH3F *hEz = (TH3F*)file->Get("hEz");

  ////////////////////////////////////////// Dx
  int Dx_xxbin = hDx->GetNbinsX();
  double Dx_xxlow = hDx->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Dx_xxhgh = hDx->GetXaxis()->GetBinUpEdge(Dx_xxbin) * 1.;
  int Dx_yybin = hDx->GetNbinsY();
  double Dx_yylow = hDx->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Dx_yyhgh = hDx->GetYaxis()->GetBinUpEdge(Dx_yybin) * 1.;
  int Dx_zzbin = hDx->GetNbinsZ();
  double Dx_zzlow = hDx->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Dx_zzhgh = hDx->GetZaxis()->GetBinUpEdge(Dx_zzbin) * 1.;
  h3_Dx = new TH3D("h3_Dx", "h3_Dx", 
		   Dx_xxbin, Dx_xxlow, Dx_xxhgh, 
		   Dx_yybin, Dx_yylow, Dx_yyhgh,
		   Dx_zzbin, Dx_zzlow, Dx_zzhgh);
  for(int i=1; i<=Dx_xxbin; i++) {
    for(int j=1; j<=Dx_yybin; j++) {
      for(int k=1; k<=Dx_zzbin; k++) {
	double content = hDx->GetBinContent(i,j,k) * 1.;
	h3_Dx->SetBinContent(i,j,k, content);
      }
    }
  }

  ////////////////////////////////////////// Dy
  int Dy_xxbin = hDy->GetNbinsX();
  double Dy_xxlow = hDy->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Dy_xxhgh = hDy->GetXaxis()->GetBinUpEdge(Dy_xxbin) * 1.;
  int Dy_yybin = hDy->GetNbinsY();
  double Dy_yylow = hDy->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Dy_yyhgh = hDy->GetYaxis()->GetBinUpEdge(Dy_yybin) * 1.;
  int Dy_zzbin = hDy->GetNbinsZ();
  double Dy_zzlow = hDy->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Dy_zzhgh = hDy->GetZaxis()->GetBinUpEdge(Dy_zzbin) * 1.;
  h3_Dy = new TH3D("h3_Dy", "h3_Dy", 
		   Dy_xxbin, Dy_xxlow, Dy_xxhgh, 
		   Dy_yybin, Dy_yylow, Dy_yyhgh,
		   Dy_zzbin, Dy_zzlow, Dy_zzhgh);
  for(int i=1; i<=Dy_xxbin; i++) {
    for(int j=1; j<=Dy_yybin; j++) {
      for(int k=1; k<=Dy_zzbin; k++) {
	double content = hDy->GetBinContent(i,j,k) * 1.;
	h3_Dy->SetBinContent(i,j,k, content);
      }
    }
  }

  ////////////////////////////////////////// Dz
  int Dz_xxbin = hDz->GetNbinsX();
  double Dz_xxlow = hDz->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Dz_xxhgh = hDz->GetXaxis()->GetBinUpEdge(Dz_xxbin) * 1.;
  int Dz_yybin = hDz->GetNbinsY();
  double Dz_yylow = hDz->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Dz_yyhgh = hDz->GetYaxis()->GetBinUpEdge(Dz_yybin) * 1.;
  int Dz_zzbin = hDz->GetNbinsZ();
  double Dz_zzlow = hDz->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Dz_zzhgh = hDz->GetZaxis()->GetBinUpEdge(Dz_zzbin) * 1.;
  h3_Dz = new TH3D("h3_Dz", "h3_Dz", 
		   Dz_xxbin, Dz_xxlow, Dz_xxhgh, 
		   Dz_yybin, Dz_yylow, Dz_yyhgh,
		   Dz_zzbin, Dz_zzlow, Dz_zzhgh);
  for(int i=1; i<=Dz_xxbin; i++) {
    for(int j=1; j<=Dz_yybin; j++) {
      for(int k=1; k<=Dz_zzbin; k++) {
	double content = hDz->GetBinContent(i,j,k) * 1.;
	h3_Dz->SetBinContent(i,j,k, content);
      }
    }
  }

  ////////////////////////////////////////// Ex
  int Ex_xxbin = hEx->GetNbinsX();
  double Ex_xxlow = hEx->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Ex_xxhgh = hEx->GetXaxis()->GetBinUpEdge(Ex_xxbin) * 1.;
  int Ex_yybin = hEx->GetNbinsY();
  double Ex_yylow = hEx->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Ex_yyhgh = hEx->GetYaxis()->GetBinUpEdge(Ex_yybin) * 1.;
  int Ex_zzbin = hEx->GetNbinsZ();
  double Ex_zzlow = hEx->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Ex_zzhgh = hEx->GetZaxis()->GetBinUpEdge(Ex_zzbin) * 1.;
  h3_Ex = new TH3D("h3_Ex", "h3_Ex", 
		   Ex_xxbin, Ex_xxlow, Ex_xxhgh, 
		   Ex_yybin, Ex_yylow, Ex_yyhgh,
		   Ex_zzbin, Ex_zzlow, Ex_zzhgh);
  for(int i=1; i<=Ex_xxbin; i++) {
    for(int j=1; j<=Ex_yybin; j++) {
      for(int k=1; k<=Ex_zzbin; k++) {
	double content = hEx->GetBinContent(i,j,k) * 1.;
	h3_Ex->SetBinContent(i,j,k, content);
      }
    }
  }

  ////////////////////////////////////////// Ey
  int Ey_xxbin = hEy->GetNbinsX();
  double Ey_xxlow = hEy->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Ey_xxhgh = hEy->GetXaxis()->GetBinUpEdge(Ey_xxbin) * 1.;
  int Ey_yybin = hEy->GetNbinsY();
  double Ey_yylow = hEy->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Ey_yyhgh = hEy->GetYaxis()->GetBinUpEdge(Ey_yybin) * 1.;
  int Ey_zzbin = hEy->GetNbinsZ();
  double Ey_zzlow = hEy->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Ey_zzhgh = hEy->GetZaxis()->GetBinUpEdge(Ey_zzbin) * 1.;
  h3_Ey = new TH3D("h3_Ey", "h3_Ey", 
		   Ey_xxbin, Ey_xxlow, Ey_xxhgh, 
		   Ey_yybin, Ey_yylow, Ey_yyhgh,
		   Ey_zzbin, Ey_zzlow, Ey_zzhgh);
  for(int i=1; i<=Ey_xxbin; i++) {
    for(int j=1; j<=Ey_yybin; j++) {
      for(int k=1; k<=Ey_zzbin; k++) {
	double content = hEy->GetBinContent(i,j,k) * 1.;
	h3_Ey->SetBinContent(i,j,k, content);
      }
    }
  }

  ////////////////////////////////////////// Ez
  int Ez_xxbin = hEz->GetNbinsX();
  double Ez_xxlow = hEz->GetXaxis()->GetBinLowEdge(1) * 1.;
  double Ez_xxhgh = hEz->GetXaxis()->GetBinUpEdge(Ez_xxbin) * 1.;
  int Ez_yybin = hEz->GetNbinsY();
  double Ez_yylow = hEz->GetYaxis()->GetBinLowEdge(1) * 1.;
  double Ez_yyhgh = hEz->GetYaxis()->GetBinUpEdge(Ez_yybin) * 1.;
  int Ez_zzbin = hEz->GetNbinsZ();
  double Ez_zzlow = hEz->GetZaxis()->GetBinLowEdge(1) * 1.;
  double Ez_zzhgh = hEz->GetZaxis()->GetBinUpEdge(Ez_zzbin) * 1.;
  h3_Ez = new TH3D("h3_Ez", "h3_Ez", 
		   Ez_xxbin, Ez_xxlow, Ez_xxhgh, 
		   Ez_yybin, Ez_yylow, Ez_yyhgh,
		   Ez_zzbin, Ez_zzlow, Ez_zzhgh);
  for(int i=1; i<=Ez_xxbin; i++) {
    for(int j=1; j<=Ez_yybin; j++) {
      for(int k=1; k<=Ez_zzbin; k++) {
	double content = hEz->GetBinContent(i,j,k) * 1.;
	h3_Ez->SetBinContent(i,j,k, content);
      }
    }
  }


  
  std::cout<<std::endl<<" Initialize Pos&Efield Correction "<<std::endl;
  std::cout<<" hDx(5,5,5) "<<h3_Dx->GetBinContent(5,5,5)<<" "<<hDx->GetBinContent(5,5,5)<<std::endl;
  std::cout<<" hDy(5,5,5) "<<h3_Dy->GetBinContent(5,5,5)<<" "<<hDy->GetBinContent(5,5,5)<<std::endl;
  std::cout<<" hDz(5,5,5) "<<h3_Dz->GetBinContent(5,5,5)<<" "<<hDz->GetBinContent(5,5,5)<<std::endl;
  std::cout<<" hEx(5,5,5) "<<h3_Ex->GetBinContent(5,5,5)<<" "<<hEx->GetBinContent(5,5,5)<<std::endl;
  std::cout<<" hEy(5,5,5) "<<h3_Ey->GetBinContent(5,5,5)<<" "<<hEy->GetBinContent(5,5,5)<<std::endl;
  std::cout<<" hEz(5,5,5) "<<h3_Ez->GetBinContent(5,5,5)<<" "<<hEz->GetBinContent(5,5,5)<<std::endl;
  std::cout<<std::endl;
  std::cout<<" maximum "<<std::endl;
  std::cout<<" hDx( maximum ) "<<h3_Dx->GetMaximum()<<" "<<hDx->GetMaximum()<<std::endl;
  std::cout<<" hDy( maximum ) "<<h3_Dy->GetMaximum()<<" "<<hDy->GetMaximum()<<std::endl;
  std::cout<<" hDz( maximum ) "<<h3_Dz->GetMaximum()<<" "<<hDz->GetMaximum()<<std::endl;
  std::cout<<" hEx( maximum ) "<<h3_Ex->GetMaximum()<<" "<<hEx->GetMaximum()<<std::endl;
  std::cout<<" hEy( maximum ) "<<h3_Ey->GetMaximum()<<" "<<hEy->GetMaximum()<<std::endl;
  std::cout<<" hEz( maximum ) "<<h3_Ez->GetMaximum()<<" "<<hEz->GetMaximum()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" minimum "<<std::endl;
  std::cout<<" hDx( minimum ) "<<h3_Dx->GetMinimum()<<" "<<hDx->GetMinimum()<<std::endl;
  std::cout<<" hDy( minimum ) "<<h3_Dy->GetMinimum()<<" "<<hDy->GetMinimum()<<std::endl;
  std::cout<<" hDz( minimum ) "<<h3_Dz->GetMinimum()<<" "<<hDz->GetMinimum()<<std::endl;
  std::cout<<" hEx( minimum ) "<<h3_Ex->GetMinimum()<<" "<<hEx->GetMinimum()<<std::endl;
  std::cout<<" hEy( minimum ) "<<h3_Ey->GetMinimum()<<" "<<hEy->GetMinimum()<<std::endl;
  std::cout<<" hEz( minimum ) "<<h3_Ez->GetMinimum()<<" "<<hEz->GetMinimum()<<std::endl;
  std::cout<<std::endl;
  std::cout<<" hDx "<<Dx_xxbin<<" "<<Dx_xxlow<<" "<<Dx_xxhgh<<", "
	   <<" hDx "<<Dx_yybin<<" "<<Dx_yylow<<" "<<Dx_yyhgh<<", "
	   <<" hDx "<<Dx_zzbin<<" "<<Dx_zzlow<<" "<<Dx_zzhgh<<std::endl;
  std::cout<<" hDy "<<Dy_xxbin<<" "<<Dy_xxlow<<" "<<Dy_xxhgh<<", "
	   <<" hDy "<<Dy_yybin<<" "<<Dy_yylow<<" "<<Dy_yyhgh<<", "
	   <<" hDy "<<Dy_zzbin<<" "<<Dy_zzlow<<" "<<Dy_zzhgh<<std::endl;
  std::cout<<" hDz "<<Dz_xxbin<<" "<<Dz_xxlow<<" "<<Dz_xxhgh<<", "
	   <<" hDz "<<Dz_yybin<<" "<<Dz_yylow<<" "<<Dz_yyhgh<<", "
	   <<" hDz "<<Dz_zzbin<<" "<<Dz_zzlow<<" "<<Dz_zzhgh<<std::endl;
  std::cout<<" hEx "<<Ex_xxbin<<" "<<Ex_xxlow<<" "<<Ex_xxhgh<<", "
	   <<" hEx "<<Ex_yybin<<" "<<Ex_yylow<<" "<<Ex_yyhgh<<", "
	   <<" hEx "<<Ex_zzbin<<" "<<Ex_zzlow<<" "<<Ex_zzhgh<<std::endl;
  std::cout<<" hEy "<<Ey_xxbin<<" "<<Ey_xxlow<<" "<<Ey_xxhgh<<", "
	   <<" hEy "<<Ey_yybin<<" "<<Ey_yylow<<" "<<Ey_yyhgh<<", "
	   <<" hEy "<<Ey_zzbin<<" "<<Ey_zzlow<<" "<<Ey_zzhgh<<std::endl;
  std::cout<<" hEz "<<Ez_xxbin<<" "<<Ez_xxlow<<" "<<Ez_xxhgh<<", "
	   <<" hEz "<<Ez_yybin<<" "<<Ez_yylow<<" "<<Ez_yyhgh<<", "
	   <<" hEz "<<Ez_zzbin<<" "<<Ez_zzlow<<" "<<Ez_zzhgh<<std::endl;

  file->Close();

  //flag_PosEfield_corr = true;

}

Point WCP::TPCParams::func_pos_SCE_correction(Point& pos){
  const double x_length = 256;
  const double y_length = 232.5;
  const double z_length = 1037;
    
  const double scale_x = 2.50/2.56;
  const double scale_y = 2.50/2.33;
  const double scale_z = 10.0/10.37;

  double p1_x = pos.x/units::cm;
  double p1_y = pos.y/units::cm;
  double p1_z = pos.z/units::cm;
  
  /////////////////////////////////////////////// p1

  p1_y += y_length/2;

  if( p1_x<0.001 )    p1_x = 0.001;
  if( p1_x>255.999 )  p1_x = 255.999 -0.001;
  if( p1_y<0.001 )    p1_y = 0.001;
  if( p1_y>232.499 )  p1_y = 232.499 -0.001;
  if( p1_z<0.001 )    p1_z = 0.001;
  if( p1_z>1036.999 ) p1_z = 1036.999 -0.001;

  double CorrWorld_p1_x = 2.50 - scale_x*(p1_x/100.0);
  double CorrWorld_p1_y = scale_y*( p1_y/100. );
  double CorrWorld_p1_z = scale_z*( p1_z/100. );

  double corr_p1_x = h3_Dx->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_y = h3_Dy->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_z = h3_Dz->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);

  p1_x = p1_x - corr_p1_x/scale_x;
  p1_y = p1_y + corr_p1_y/scale_y - y_length/2;
  p1_z = p1_z + corr_p1_z/scale_z;

  Point p(p1_x*units::cm, p1_y*units::cm, p1_z*units::cm);
  return p;
}

double WCP::TPCParams::func_dx_after_Pos_Efield_SCE_correction(double p1_x, double p1_y, double p1_z, double pA_x, double pA_y, double pA_z, double p2_x, double p2_y, double p2_z)// unit:: cm,  p1 --> pA --> p2
{
  double result = 1;

  const double x_length = 256;
  const double y_length = 232.5;
  const double z_length = 1037;
    
  const double scale_x = 2.50/2.56;
  const double scale_y = 2.50/2.33;
  const double scale_z = 10.0/10.37;

  /////////////////////////////////////////////// p1

  p1_y += y_length/2;

  if( p1_x<0.001 )    p1_x = 0.001;
  if( p1_x>255.999 )  p1_x = 255.999 -0.001;
  if( p1_y<0.001 )    p1_y = 0.001;
  if( p1_y>232.499 )  p1_y = 232.499 -0.001;
  if( p1_z<0.001 )    p1_z = 0.001;
  if( p1_z>1036.999 ) p1_z = 1036.999 -0.001;

  double CorrWorld_p1_x = 2.50 - scale_x*(p1_x/100.0);
  double CorrWorld_p1_y = scale_y*( p1_y/100. );
  double CorrWorld_p1_z = scale_z*( p1_z/100. );

  double corr_p1_x = h3_Dx->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_y = h3_Dy->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_z = h3_Dz->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);

  p1_x = p1_x - corr_p1_x/scale_x;
  p1_y = p1_y + corr_p1_y/scale_y;
  p1_z = p1_z + corr_p1_z/scale_z;

  /////////////////////////////////////////////// pA

  pA_y += y_length/2;

  if( pA_x<0.001 )    pA_x = 0.001;
  if( pA_x>255.999 )  pA_x = 255.999 -0.001;
  if( pA_y<0.001 )    pA_y = 0.001;
  if( pA_y>232.499 )  pA_y = 232.499 -0.001;
  if( pA_z<0.001 )    pA_z = 0.001;
  if( pA_z>1036.999 ) pA_z = 1036.999 -0.001;

  double CorrWorld_pA_x = 2.50 - scale_x*(pA_x/100.0);
  double CorrWorld_pA_y = scale_y*( pA_y/100. );
  double CorrWorld_pA_z = scale_z*( pA_z/100. );

  double corr_pA_x = h3_Dx->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_y = h3_Dy->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_z = h3_Dz->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);

  pA_x = pA_x - corr_pA_x/scale_x;
  pA_y = pA_y + corr_pA_y/scale_y;
  pA_z = pA_z + corr_pA_z/scale_z;

  /////////////////////////////////////////////// p2

  p2_y += y_length/2;

  if( p2_x<0.001 )    p2_x = 0.001;
  if( p2_x>255.999 )  p2_x = 255.999 -0.001;
  if( p2_y<0.001 )    p2_y = 0.001;
  if( p2_y>232.499 )  p2_y = 232.499 -0.001;
  if( p2_z<0.001 )    p2_z = 0.001;
  if( p2_z>1036.999 ) p2_z = 1036.999 -0.001;

  double CorrWorld_p2_x = 2.50 - scale_x*(p2_x/100.0);
  double CorrWorld_p2_y = scale_y*( p2_y/100. );
  double CorrWorld_p2_z = scale_z*( p2_z/100. );

  double corr_p2_x = h3_Dx->Interpolate(CorrWorld_p2_x, CorrWorld_p2_y, CorrWorld_p2_z);
  double corr_p2_y = h3_Dy->Interpolate(CorrWorld_p2_x, CorrWorld_p2_y, CorrWorld_p2_z);
  double corr_p2_z = h3_Dz->Interpolate(CorrWorld_p2_x, CorrWorld_p2_y, CorrWorld_p2_z);

  p2_x = p2_x - corr_p2_x/scale_x;
  p2_y = p2_y + corr_p2_y/scale_y;
  p2_z = p2_z + corr_p2_z/scale_z;

  ////////////////////////////////////////////////

  double result_p1 = sqrt( pow(p1_x-pA_x,2) + pow(p1_y-pA_y,2) + pow(p1_z-pA_z,2) );
  double result_p2 = sqrt( pow(p2_x-pA_x,2) + pow(p2_y-pA_y,2) + pow(p2_z-pA_z,2) );
  result = ( result_p1 + result_p2 )/2;

  return result;
}




double WCP::TPCParams::func_dx_after_Pos_Efield_SCE_correction(double p1_x, double p1_y, double p1_z, double pA_x, double pA_y, double pA_z)// unit:: cm, p1 ---> pA
{
  double result = 1;

  const double x_length = 256;
  const double y_length = 232.5;
  const double z_length = 1037;
    
  const double scale_x = 2.50/2.56;
  const double scale_y = 2.50/2.33;
  const double scale_z = 10.0/10.37;

  /////////////////////////////////////////////// p1

  p1_y += y_length/2;

  if( p1_x<0.001 )    p1_x = 0.001;
  if( p1_x>255.999 )  p1_x = 255.999 -0.001;
  if( p1_y<0.001 )    p1_y = 0.001;
  if( p1_y>232.499 )  p1_y = 232.499 -0.001;
  if( p1_z<0.001 )    p1_z = 0.001;
  if( p1_z>1036.999 ) p1_z = 1036.999 -0.001;

  double CorrWorld_p1_x = 2.50 - scale_x*(p1_x/100.0);
  double CorrWorld_p1_y = scale_y*( p1_y/100. );
  double CorrWorld_p1_z = scale_z*( p1_z/100. );

  double corr_p1_x = h3_Dx->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_y = h3_Dy->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);
  double corr_p1_z = h3_Dz->Interpolate(CorrWorld_p1_x, CorrWorld_p1_y, CorrWorld_p1_z);

  p1_x = p1_x - corr_p1_x/scale_x;
  p1_y = p1_y + corr_p1_y/scale_y;
  p1_z = p1_z + corr_p1_z/scale_z;

  /////////////////////////////////////////////// pA

  pA_y += y_length/2;

  if( pA_x<0.001 )    pA_x = 0.001;
  if( pA_x>255.999 )  pA_x = 255.999 -0.001;
  if( pA_y<0.001 )    pA_y = 0.001;
  if( pA_y>232.499 )  pA_y = 232.499 -0.001;
  if( pA_z<0.001 )    pA_z = 0.001;
  if( pA_z>1036.999 ) pA_z = 1036.999 -0.001;

  double CorrWorld_pA_x = 2.50 - scale_x*(pA_x/100.0);
  double CorrWorld_pA_y = scale_y*( pA_y/100. );
  double CorrWorld_pA_z = scale_z*( pA_z/100. );

  double corr_pA_x = h3_Dx->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_y = h3_Dy->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_z = h3_Dz->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);

  pA_x = pA_x - corr_pA_x/scale_x;
  pA_y = pA_y + corr_pA_y/scale_y;
  pA_z = pA_z + corr_pA_z/scale_z;

  ////////////////////////////////////////////////

  double result_p1 = sqrt( pow(p1_x-pA_x,2) + pow(p1_y-pA_y,2) + pow(p1_z-pA_z,2) );
  result = result_p1;

  return result;
}



double WCP::TPCParams::func_dQdx_after_Pos_Efield_SCE_correction(double pA_x, double pA_y, double pA_z, double dQ, double dx)
{
  double result = 1;

  const double x_length = 256;
  const double y_length = 232.5;
  const double z_length = 1037;
    
  const double scale_x = 2.50/2.56;
  const double scale_y = 2.50/2.33;
  const double scale_z = 10.0/10.37;

  const double E0 = 0.2739;  
  const double alpha_ArgoNeut = 0.93;
  const double beta_ArgoNeut  = 0.212;

  /////////////////////////////////////////////// pA

  pA_y += y_length/2;

  if( pA_x<0.001 )    pA_x = 0.001;
  if( pA_x>255.999 )  pA_x = 255.999 -0.001;
  if( pA_y<0.001 )    pA_y = 0.001;
  if( pA_y>232.499 )  pA_y = 232.499 -0.001;
  if( pA_z<0.001 )    pA_z = 0.001;
  if( pA_z>1036.999 ) pA_z = 1036.999 -0.001;

  double CorrWorld_pA_x = 2.50 - scale_x*(pA_x/100.0);
  double CorrWorld_pA_y = scale_y*( pA_y/100. );
  double CorrWorld_pA_z = scale_z*( pA_z/100. );

  double corr_pA_x = h3_Ex->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_y = h3_Ey->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);
  double corr_pA_z = h3_Ez->Interpolate(CorrWorld_pA_x, CorrWorld_pA_y, CorrWorld_pA_z);

  double Ex = E0 - corr_pA_x;
  double Ey = corr_pA_y;
  double Ez = corr_pA_z;      
  double Efield = sqrt( Ex*Ex + Ey*Ey + Ez*Ez );

  double dQdx_before = dQ/dx;
  
  double val_dEdx_from_dQdx = func_dEdx_from_dQdx( dQdx_before, Efield, alpha_ArgoNeut, beta_ArgoNeut );
  double val_dQdx_from_dEdx = func_dQdx_from_dEdx_by_ArgoNeut_model( val_dEdx_from_dQdx, E0, alpha_ArgoNeut, beta_ArgoNeut );
  result = val_dQdx_from_dEdx;

  return result;
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
