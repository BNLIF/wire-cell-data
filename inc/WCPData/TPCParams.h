#ifndef WIRECELL_TPCPARAMS_H
#define WIRECELL_TPCPARAMS_H

#include "TGraph.h"
#include "TString.h"
#include "WCPData/Point.h"

#include "TH3.h"

namespace WCP{
  class TPCParams {
    double m_pitch_u; // wire pitch u
    double m_pitch_v; // wire pitch v
    double m_pitch_w; // wire pitch w
    
    double m_ts_width; // time slice width 2 us * 1.6 mm/us ~ 3.2 mm

    double m_angle_u;
    double m_angle_v;
    double m_angle_w;

    double first_u_dis;
    double first_v_dis;
    double first_w_dis;

    int nrebin;
    int time_offset;

    bool flag_corr;
    TGraph *gu, *gv, *gw;

    TGraph *g_proton, *g_muon, *g_pion, *g_kaon, *g_electron;

    double mass_proton;
    double mass_neutron;
    double mass_muon;
    double mass_pion;
    double mass_neutral_pion;
    double mass_kaon;
    double mass_electron;
    
    TGraph *g_proton_r2ke, *g_muon_r2ke, *g_pion_r2ke, *g_kaon_r2ke, *g_electron_r2ke;

    
    double electron_lifetime; // value in ms ...

    
    bool flag_PosEfield_corr;
    /* TH3F *hDx; */
    /* TH3F *hDy; */
    /* TH3F *hDz; */
    /* TH3F *hEx; */
    /* TH3F *hEy; */
    /* TH3F *hEz; */    
    TH3D *h3_Dx;
    TH3D *h3_Dy;
    TH3D *h3_Dz;
    TH3D *h3_Ex;
    TH3D *h3_Ey;
    TH3D *h3_Ez;
  public:
    TPCParams();
    ~TPCParams();
    
    double get_mass_proton(){return mass_proton;};
    double get_mass_neutron(){return mass_neutron;};
    double get_mass_kaon(){return mass_kaon;};
    double get_mass_pion(){return mass_pion;};
    double get_mass_muon(){return mass_muon;};
    double get_mass_pi0(){return mass_neutral_pion;};
    double get_mass_electron(){return mass_electron;};
    
    // Position and E-field correction for SCE
    void init_Pos_Efield_SCE_correction(TString filename="input_data_files/SCEoffsets_dataDriven_combined_bkwd_Jan18.root");
    bool get_flag_PosEfield_corr() { return flag_PosEfield_corr; }
    double func_dx_after_Pos_Efield_SCE_correction(double p1_x, double p1_y, double p1_z, double pA_x, double pA_y, double pA_z, double p2_x, double p2_y, double p2_z);// unit:: cm,  p1 --> pA --> p2
    double func_dx_after_Pos_Efield_SCE_correction(double p1_x, double p1_y, double p1_z, double pA_x, double pA_y, double pA_z);// unit:: cm,  p1 --> pA

    void set_flag_PosEfield_corr(bool val){flag_PosEfield_corr = val;};
									
    Point func_pos_SCE_correction(Point& pos);
    

    double func_dQdx_from_dEdx_by_ArgoNeut_model(double dEdx, double e, double alpha, double beta)
{
  double result = log(alpha + beta/1.38/e*dEdx)/(23.6e-6*beta/1.38/e);
  return result;
}

    double func_dEdx_from_dQdx(double dQdx, double e, double alpha, double beta)
{
  double result = exp( dQdx * beta * 23.6e-6/1.38/e ) - alpha;
  result = result/( beta/1.38/e );
  return result;
}

    double func_dQdx_after_Pos_Efield_SCE_correction(double pA_x, double pA_y, double pA_z, double dQ, double dx);

    //set/get u first dis
    void set_first_u_dis(double p){ first_u_dis = p;}
    double get_first_u_dis(){return first_u_dis;}

    //set/get u first dis
    void set_first_v_dis(double p){ first_v_dis = p;}
    double get_first_v_dis(){return first_v_dis;}

    //set/get u first dis
    void set_first_w_dis(double p){ first_w_dis = p;}
    double get_first_w_dis(){return first_w_dis;}
    
    // set/get u pitches
    void set_pitch_u(double p) { m_pitch_u = p; }
    double get_pitch_u() { return m_pitch_u; }

    // set/get v pitches
    void set_pitch_v(double p) { m_pitch_v = p; }
    double get_pitch_v() { return m_pitch_v; }
    
    // set/get w pitches
    void set_pitch_w(double p) { m_pitch_w = p; }
    double get_pitch_w() { return m_pitch_w; }
    
    double get_pitch(){return (m_pitch_u + m_pitch_v + m_pitch_w)/3.;};
    
    void set_ts_width(double p){ m_ts_width = p;}
    double get_ts_width(){return m_ts_width;}

    void set_angle_u(double p){m_angle_u = p;}
    double get_angle_u(){return m_angle_u;};

    void set_angle_v(double p){m_angle_v = p;}
    double get_angle_v(){return m_angle_v;};

    void set_angle_w(double p){m_angle_w = p;}
    double get_angle_w(){return m_angle_w;};
    

    void set_time_offset(int p){time_offset = p;};
    int get_time_offset(){return time_offset;};

    void set_nrebin(int p){nrebin = p;};
    int get_nrebin(){return nrebin;};


    void set_electron_lifetime(double p){electron_lifetime = p;};
    double get_electron_lifetime(){return electron_lifetime;};
    double get_attenuation_ratio(double drift_time); // in ms ...

    
    bool get_flag_corr(){return flag_corr;};
    
    // etc for other parameters you need
    void init_corr_files(TString file_u="input_data_files/calib_u_corr.txt", int ndata_u = 2401, TString file_v="input_data_files/calib_v_corr.txt", int ndata_v = 2401, TString file_w="input_data_files/calib_w_corr.txt", int ndata_w = 3457);
    double get_corr_factor(WCP::Point& p, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw);

    void init_PID_dq_dx(TString filename = "input_data_files/stopping_ave_dQ_dx.root", TString filename1 = "input_data_files/ave_range_to_kenergy.root");
    
    TGraph* get_pion_dq_dx(){return g_pion;};
    TGraph* get_proton_dq_dx(){return g_proton;};
    TGraph* get_muon_dq_dx(){return g_muon;};
    TGraph* get_kaon_dq_dx(){return g_kaon;};
    TGraph *get_electron_dq_dx(){return g_electron;};

    TGraph* get_pion_r2ke(){return g_pion_r2ke;};
    TGraph* get_proton_r2ke(){return g_proton_r2ke;};
    TGraph* get_muon_r2ke(){return g_muon_r2ke;};
    TGraph* get_kaon_r2ke(){return g_kaon_r2ke;};
    TGraph *get_electron_r2ke(){return g_electron_r2ke;};
  };
 }


#endif 
