
double PR3DCluster::cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, double w_center, double w_sigma, int flag, double nsigma){
  // flag  = 0 no boundary effect, pure Gaussian, time or collection plane
  // flag  = 1 taking into account boundary effect for induction plane
  double result = 0;
  // time domain ...
  if (fabs((tbin + 0.5)-t_center) <= nsigma * t_sigma &&
      fabs((wbin + 0.5)-w_center) <= nsigma * w_sigma){
    // time dimension ...
    result = 0.5*(std::erf((tbin+1-t_center)/sqrt(2.)/t_sigma)-std::erf((tbin-t_center)/sqrt(2.)/t_sigma));

    if (flag ==0){
      result *= 0.5*(std::erf((wbin+1-w_center)/sqrt(2.)/w_sigma)-std::erf((wbin-w_center)/sqrt(2.)/w_sigma));
    }else if (flag==1){

      double x2 = wbin + 1.5;
      double x1 = wbin + 0.5;
      double x0 = w_center;
      double content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);
      //w1 = 0.5;
      // std::cout << w1 << std::endl;
      
      /* double w1 = (1./2./sqrt(2*3.1415926)/pow(x1-x2,2) * */
      /* 	    (2*exp(-pow(x0-x1,2)/2./pow(w_sigma,2)) * w_sigma *(x0+x1-2*x2) - 2 * exp(-pow(x0-x2,2)/2./pow(w_sigma,2))*w_sigma*(x0-x2)- */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x1-x0)/sqrt(2.)/w_sigma) + */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x2-x0)/sqrt(2.)/w_sigma)))/content1; */

      x2 = wbin + 0.5;
      x1 = wbin - 0.5;
      double content2 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w2 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);
      // std::cout << w2 << std::endl
      /* double w2 = (1./2./sqrt(2*3.1415926)/pow(x1-x2,2) * */
      /* 	    (2*exp(-pow(x0-x1,2)/2./pow(w_sigma,2)) * w_sigma *(x0+x1-2*x2) - 2 * exp(-pow(x0-x2,2)/2./pow(w_sigma,2))*w_sigma*(x0-x2)- */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x1-x0)/sqrt(2.)/w_sigma) + */
      /* 	     sqrt(2*3.1415926)*(pow(w_sigma,2)+pow(x0-x2,2))*std::erf((x2-x0)/sqrt(2.)/w_sigma)))/content2; */

      /* if (w1 > 0.5) { */
      /* 	w1 = sqrt(w1); */
      /* }else if (w1<0.5){ */
      /* 	w1 = 1-sqrt(1-w1); */
      /* } */
      /*  if (w2 > 0.5) { */
      /* 	w2 = sqrt(w1); */
      /*  }else if (w2 < 0.5){ */
      /* 	w2 = 1-sqrt(1-w2); */
      /* } */
      
      result *= (content1 * w1+content2 * (1-w2));
    }else if (flag==2){
      double sum = 0;
      
      double x2 = wbin + 1.0;
      double x1 = wbin + 0.5;
      double x0 = w_center;
      double content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      double w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum += content1 * (0.545 + 0.697*w1);

      x2 = wbin + 1.5;
      x1 = wbin + 1.0;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.11364 + 0.1 * w1);

      x2 = wbin + 0.5;
      x1 = wbin;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.545 + 0.697*(1-w1));

      x2 = wbin;
      x1 = wbin - 0.5;
      content1 = 0.5*(std::erf((x2-x0)/sqrt(2.)/w_sigma)-std::erf((x1-x0)/sqrt(2.)/w_sigma));
      w1 = - pow(w_sigma,2)/(-1)/sqrt(2.*3.1415926)/w_sigma
      	*(exp(-pow(x0-x2,2)/2./pow(w_sigma,2))-exp(-pow(x0-x1,2)/2./pow(w_sigma,2)))
      	/(0.5*std::erf((x2-x0)/sqrt(2.)/w_sigma) - 0.5*std::erf((x1-x0)/sqrt(2.)/w_sigma))
      	+ (x0-x2)/(-1);

      sum+= content1 *(0.11364 + 0.1 * (1-w1));
      
      result *= sum;
    }
  }
  return result;
}


void PR3DCluster::dQ_dx_fit(std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, double flash_time){

  // Need to take into account the time, so one can properly calculate X value for diffusion ...
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double time_slice_width = mp.get_ts_width();
  int time_offset = mp.get_time_offset();
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
  flash_time = flash_time - time_offset*units::microsecond; // us
  // given an x position value, we want to convert this to a drift time 
  // pos_x/time_slice_width * nrebin * 0.5*units::us // us
  // difference between these two numbers are the time in us ... 
  
  // Now figure out the diffusion coefficients
  // these are current  numbers in WCT, not sure what would be the values for data ... 
  double DL = 6.4 * pow(units::cm,2)/units::second ;
  double DT = 9.8 * pow(units::cm,2)/units::second ;

  // these are the transverse broading due to software filters in the wire dimension
  // these should be quadrature added to the 
  double col_sigma_w_T = 0.188060 * pitch_w*0.2; // units::mm
  double ind_sigma_u_T = 0.402993 * pitch_u*0.2; // units::mm
  double ind_sigma_v_T = 0.402993 * pitch_v*0.2; // units::mm
  
  // this is the longitudinal filters in the time dimension ...
  double add_sigma_L = 1.428249  * time_slice_width / nrebin / 0.5; // units::mm 


  // Now start the fit ... 
  Point reco_pos(150*units::cm,30*units::cm,600*units::cm);
  reco_pos.x = (reco_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ;
  //  std::cout << reco_pos.x/units::cm << " " << reco_pos.y/units::cm << " " << reco_pos.z/units::cm << std::endl;

  
  // let's convert these into time and wire numbers 

  // T U, V, W, wires ...
  double central_T = offset_t + slope_xt * reco_pos.x ;
  double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
  double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400;
  double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800;
  
  // start to work out the diffusion coefficients ...
  double drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ;

  //  std::cout << drift_time/units::microsecond << std::endl;
  double diff_sigma_L = sqrt(2* DL * drift_time);
  double diff_sigma_T = sqrt(2* DT * drift_time);
  
  double sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width;
  double sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u;
  double sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v;
  double sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w;
  //std::cout << diff_sigma_T/units::cm << " " << sigma_T_u/units::cm << std::endl;
  /* std::cout << drift_time/units::microsecond << " " << diff_sigma_L/units::cm << " " << diff_sigma_L/time_slice_width << " " << add_sigma_L/time_slice_width << std::endl; */
  /* std::cout << diff_sigma_T/units::mm << " " << col_sigma_w_T/units::mm << std::endl; */
  //  std::cout << central_T << " " << central_W << " " << central_U << " " << central_V << std::endl; 
  /* sigma_T_w = 0.5785; */
  /* sigma_L = 0.9818; */
  // Now, need to add the calculations ... 
 
  std::vector<int> proj_channel;
  std::vector<int> proj_timeslice;
  std::vector<int> proj_charge;
  get_projection(proj_channel,proj_timeslice,proj_charge, global_wc_map);

  // condense the information ...
  std::map<std::pair<int,int>, double> proj_data_map;
  for (size_t i=0;i!=proj_charge.size(); i++){
    auto it = proj_data_map.find(std::make_pair(proj_channel.at(i),proj_timeslice.at(i)));
    if (it == proj_data_map.end()){
      proj_data_map[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = proj_charge.at(i);
    }else{
      it->second += proj_charge.at(i);
    }
  }
  
  for (auto it = proj_data_map.begin(); it!=proj_data_map.end(); it++){
    if (it->first.first>=4800){
      double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_W, sigma_T_w,0,4)*500000;
      std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; 
    }else if (it->first.first>=2400){
      double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_V, sigma_T_v,0,4)*500000;
      std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl;
    }else{
      double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_U, sigma_T_u,0,4)*500000;
      std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl;
    }
  }
 
  // these should be the expected values:
  // U, V, W, T,  1252.01, 3819.54, 6799.62, 1485.81
  // reco position 1252.02, 3819.63, 6799.67, 1485.79
  // good ...

  
  //  
  //  std::cout << first_t_dis/units::cm << " " << time_slice_width/units::cm << std::endl; 
  // double offset_x = (flash_time - time_offset)*2./nrebin*time_slice_width;

  
  
  //std::cout << time_slice_width/units::cm << " " << offset_x << std::endl;
  


  // input: fine tracking trajectory ... 
  // See the collected 2D projection points ...
  
  // Need to take into account the software filters regarding the induction and collection ...

  
}
