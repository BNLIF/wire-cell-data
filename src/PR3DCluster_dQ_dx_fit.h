
void PR3DCluster::dQ_dx_fit(double flash_time){

  // Need to take into account the time, so one can properly calculate X value for diffusion ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_tick_offset = mp.get_time_tick_offset();
  int nrebin = mp.get_nrebin();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w + 0.5;
  double offset_u = -first_u_dis/pitch_u + 0.5;
  double offset_v = -first_v_dis/pitch_v + 0.5;

  double first_t_dis = path_wcps.front().mcell->GetTimeSlice()*time_slice_width - path_wcps.front().x;
  double slope_xt = 1./time_slice_width;
  double offset_t =  first_t_dis/time_slice_width + 0.5;
  
  // get the correct flash time matching TPC side
  flash_time = flash_time - time_tick_offset*0.5*units::microsecond; // us
  // given an x position value, we want to convert this to a drift time 
  // pos_x/time_slice_width * nrebin * 0.5*units::us // us
  // difference between these two numbers are the time in us ... 
  
  // Now figure out the diffusion coefficients
  // these are current  numbers in WCT, not sure what would be the values for data ... 
  double DL = 7.2 * pow(units::cm,2)/units::second ;
  double DT = 12.0 * pow(units::cm,2)/units::second ;

  // these are the transverse broading due to software filters in the wire dimension
  // these should be quadrature added to the 
  double col_sigma_w_T = 0.188060 * pitch_w; // units::mm
  double ind_sigma_u_T = 0.402993 * pitch_u; // units::mm
  double ind_sigma_v_T = 0.402993 * pitch_v; // units::mm
  
  // these is the longitudinal filters in the time dimension ...
  double sigma_L = 1.13656 * time_slice_width * units::cm / nrebin / 0.5; // units::mm 


  // Now start the fit ... 
  Point reco_pos(150*units::cm,30*units::cm,600*units::cm);
  reco_pos.x = (reco_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ;
  //  std::cout << reco_pos.x/units::cm << " " << reco_pos.y/units::cm << " " << reco_pos.z/units::cm << std::endl;

  
  // let's convert these into time and wire numbers 

  // T U, V, W, wires ...
 double central_T = offset_t + slope_xt * reco_pos.x ;
 double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
 double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z);
 double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z);
 
 
  // these should be the expected values:
  // U, V, W, T,  1252.01, 3819.54, 6799.62, 1485.81
  // reco position 1252.02, 3819.63, 6799.67, 1485.79
  // good ...

  
  //  
  //  std::cout << first_t_dis/units::cm << " " << time_slice_width/units::cm << std::endl; 
  // double offset_x = (flash_time - time_tick_offset)*2./nrebin*time_slice_width;

  
  
  //std::cout << time_slice_width/units::cm << " " << offset_x << std::endl;
  


  // input: fine tracking trajectory ... 
  // See the collected 2D projection points ...
  
  // Need to take into account the software filters regarding the induction and collection ...

  
}
