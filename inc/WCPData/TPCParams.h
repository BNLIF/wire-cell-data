#ifndef WIRECELL_TPCPARAMS_H
#define WIRECELL_TPCPARAMS_H

#include "TGraph.h"
#include "TString.h"
#include "WCPData/Point.h"

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
    
  public:
    TPCParams();
    ~TPCParams();

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

    bool get_flag_corr(){return flag_corr;};
    
    // etc for other parameters you need
    void init_corr_files(TString file_u="input_data_files/calib_u_corr.txt", int ndata_u = 2401, TString file_v="input_data_files/calib_v_corr.txt", int ndata_v = 2401, TString file_w="input_data_files/calib_w_corr.txt", int ndata_w = 3457);
    double get_corr_factor(WCP::Point& p, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw);
  };
 }


#endif 
