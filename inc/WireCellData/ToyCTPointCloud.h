#ifndef WireCell_TOYCTPOINTCLOUD_H
#define WireCell_TOYCTPOINTCLOUD_H

#include "WireCellNanoflann/nanoflann.h"
#include "WireCellData/WCPointCloud.h"

#include <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
nanoflann::L2_Simple_Adaptor<double, WireCell::CTPointCloud<double> > ,
  WireCell::CTPointCloud<double>,
  2 /* dim */
  > my_kd_tree_ct_t;

namespace WireCell{
  class ToyCTPointCloud {
  public:
    ToyCTPointCloud(int u_min_ch, int u_max_ch, int v_min_ch, int v_max_ch, int w_min_ch, int w_max_ch, double offset_t, double offset_u, double offset_v, double offset_w, double slope_t, double slope_u, double slope_v, double slope_w, double angle_u = 1.0472, double angle_v = -1.0472, double angle_w = 0);
    ~ToyCTPointCloud();

    void AddPoint(int ch, int time_sice, int charge, int charge_err);
    void AddPoints(std::vector<int> *timesliceId, std::vector<std::vector<int>> *timesliceChannel, std::vector<std::vector<int>> *raw_charge , std::vector<std::vector<int>> *raw_charge_err);
    void AddDeadChs(std::map<int,std::pair<double,double> >& dead_u_index, std::map<int,std::pair<double,double> >& dead_v_index, std::map<int,std::pair<double,double> >& dead_w_index);
    void build_kdtree_index();

    int get_num_points(int plane);
    WireCell::CTPointCloud<double>& get_cloud(int plane);

    bool is_good_point(WireCell::Point& p, double radius = 0.6*units::cm, int ch_range = 1, int allowed_bad = 1);
    
    WireCell::CTPointCloud<double> get_closest_points(WireCell::Point& p, double radius, int plane);
    bool get_closest_dead_chs(WireCell::Point& p, int plane, int ch_range=1);

    void Print(WireCell::Point &p);

    std::vector<int> convert_3Dpoint_time_ch(WireCell::Point& p);

    std::pair<int,int> get_uch_limits(){return std::make_pair(u_min_ch, u_max_ch);};
    std::pair<int,int> get_vch_limits(){return std::make_pair(v_min_ch, v_max_ch);};
    std::pair<int,int> get_wch_limits(){return std::make_pair(w_min_ch, w_max_ch);};
    
  protected:

    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, double radius, int plane);

    std::map<int, std::pair<double, double> > dead_uchs;
    std::map<int, std::pair<double, double> > dead_vchs;
    std::map<int, std::pair<double, double> > dead_wchs;
    
    double angle_u, angle_v, angle_w; // wire angles 
    int u_min_ch, u_max_ch; // channel range for U
    int v_min_ch, v_max_ch; // channel range for V
    int w_min_ch, w_max_ch; // channel range for W 
    
    // convert time into a position
    // (time_slice - offset_t) / slope_t = position_x 
    double offset_t, slope_t;
    
    // convert u wire number into a position
    // (u_index -offset_u) / slope_u = position_u
    double offset_u, slope_u;
    // convert v wire number into a position
    double offset_v, slope_v;
    // convert w wire number into a position 
    double offset_w, slope_w;
      
    WireCell::CTPointCloud<double> cloud_u;
    my_kd_tree_ct_t *index_u;
    WireCell::CTPointCloud<double> cloud_v;
    my_kd_tree_ct_t *index_v;
    WireCell::CTPointCloud<double> cloud_w;
    my_kd_tree_ct_t *index_w;
  };
}


#endif

