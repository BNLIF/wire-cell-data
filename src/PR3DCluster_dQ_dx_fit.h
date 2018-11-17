
double PR3DCluster::cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag, double nsigma){
  double result = 0;
  double result1 = 0;

  for (size_t i=0;i!=t_centers.size();i++){
    result += cal_gaus_integral(tbin,wbin,t_centers.at(i), t_sigmas.at(i), w_centers.at(i), w_sigmas.at(i),flag,nsigma) * weights.at(i);
    result1 += weights.at(i);
  }

  result /= result1;
  
  return result;
}

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
  // there is some cancellation of this effect with the effect of using average field response in the deconvolution ... 0-
  double col_sigma_w_T = 0.188060 * pitch_w*0.2; // units::mm
  double ind_sigma_u_T = 0.402993 * pitch_u*0.3; // units::mm
  double ind_sigma_v_T = 0.402993 * pitch_v*0.5; // units::mm
  
  // this is the longitudinal filters in the time dimension ...
  double add_sigma_L = 1.428249  * time_slice_width / nrebin / 0.5; // units::mm 

  // Point-like case ... 
  // these should be the expected values:
  // U, V, W, T,  1252.01, 3819.54, 6799.62, 1485.81
  // reco position 1252.02, 3819.63, 6799.67, 1485.79
  // good ...
  /* // Now start the fit ...  */
  /* Point reco_pos(150*units::cm-3*units::mm/2./sqrt(3.),30*units::cm-3*units::mm/2./sqrt(3.),600*units::cm-3*units::mm/2./sqrt(3.)); */
  /* reco_pos.x = (reco_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
  /* //  std::cout << reco_pos.x/units::cm << " " << reco_pos.y/units::cm << " " << reco_pos.z/units::cm << std::endl; */
  /* // let's convert these into time and wire numbers  */

  /* // T U, V, W, wires ... */
  /* double central_T = offset_t + slope_xt * reco_pos.x ; */
  /* double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z); */
  /* double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400; */
  /* double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800; */
  
  /* // start to work out the diffusion coefficients ... */
  /* double drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ; */
  /* //  std::cout << drift_time/units::microsecond << std::endl; */
  /* double diff_sigma_L = sqrt(2* DL * drift_time); */
  /* double diff_sigma_T = sqrt(2* DT * drift_time); */
  
  /* double sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width; */
  /* double sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u; */
  /* double sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v; */
  /* double sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w; */
  /* //std::cout << diff_sigma_T/units::cm << " " << sigma_T_u/units::cm << std::endl; */
  /* /\* std::cout << drift_time/units::microsecond << " " << diff_sigma_L/units::cm << " " << diff_sigma_L/time_slice_width << " " << add_sigma_L/time_slice_width << std::endl; *\/ */
  /* /\* std::cout << diff_sigma_T/units::mm << " " << col_sigma_w_T/units::mm << std::endl; *\/ */
  /* std::cout << central_T << " " << central_W << " " << central_U << " " << central_V << std::endl;  */



  /* // short track segmentation ...  */
  /* // Now, a segment ... */
  /* std::vector<double> centers_U ; */
  /* std::vector<double> centers_V ; */
  /* std::vector<double> centers_W ; */
  /* std::vector<double> centers_T ; */
  /* std::vector<double> sigmas_T; */
  /* std::vector<double> sigmas_U; */
  /* std::vector<double> sigmas_V; */
  /* std::vector<double> sigmas_W; */
  
  /* for (int i=0;i!=10;i++){ */
  /*   Point reco_pos(150*units::cm-2.9*units::mm/sqrt(3.)/9.*i,30*units::cm-2.9*units::mm/sqrt(3.)/9.*i,600*units::cm-2.9*units::mm/sqrt(3.)/9.*i); */
  /*   reco_pos.x = (reco_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
  /*   // T U, V, W, wires ... */
  /*   double central_T = offset_t + slope_xt * reco_pos.x ; */
  /*   double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z); */
  /*   double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400; */
  /*   double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800; */

  /*   // start to work out the diffusion coefficients ... */
  /*   double drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ; */
  /*   //  std::cout << drift_time/units::microsecond << std::endl; */
  /*   double diff_sigma_L = sqrt(2* DL * drift_time); */
  /*   double diff_sigma_T = sqrt(2* DT * drift_time); */
    
  /*   double sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width; */
  /*   double sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u; */
  /*   double sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v; */
  /*   double sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w; */

  /*   centers_U.push_back(central_U); */
  /*   centers_V.push_back(central_V); */
  /*   centers_W.push_back(central_W); */
  /*   centers_T.push_back(central_T); */

  /*   sigmas_U.push_back(sigma_T_u); */
  /*   sigmas_V.push_back(sigma_T_v); */
  /*   sigmas_W.push_back(sigma_T_w); */
  /*   sigmas_T.push_back(sigma_L); */
  /* } */


  
  


  
  /* for (auto it = proj_data_map.begin(); it!=proj_data_map.end(); it++){ */
  /*   /\* if (it->first.first>=4800){ *\/ */
  /*   /\*   double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_W, sigma_T_w,0,4)*500000; *\/ */
  /*   /\*   std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl;  *\/ */
  /*   /\* }else if (it->first.first>=2400){ *\/ */
  /*   /\*   double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_V, sigma_T_v,0,4)*500000; *\/ */
  /*   /\*   std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; *\/ */
  /*   /\* }else{ *\/ */
  /*   /\*   double value = cal_gaus_integral(it->first.second, it->first.first,central_T, sigma_L, central_U, sigma_T_u,0,4)*500000; *\/ */
  /*   /\*   std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; *\/ */
  /*   /\* } *\/ */

  /*   if (it->first.first>=4800){ */
  /*     double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_W, sigmas_W,0,4)*5000*30; */
  /*     std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; */
  /*   }else if (it->first.first>=2400){ */
  /*     double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_V, sigmas_V,0,4)*5000*30; */
  /*     std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; */
  /*   }else{ */
  /*     double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_U, sigmas_U,0,4)*5000*30; */
  /*     std::cout << it->first.first << " " << it->first.second << " " << it->second << " " << value <<  std::endl; */
  /*   } */
  /* } */
 
  // Now, need to add the calculations ...
  std::vector<int> proj_channel;
  std::vector<int> proj_timeslice;
  std::vector<int> proj_charge;
  std::vector<int> proj_charge_err;
  get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, global_wc_map);
  // condense the information ...
  std::map<std::pair<int,int>, std::pair<double,double> > proj_data_u_map, proj_data_v_map, proj_data_w_map;
  for (size_t i=0;i!=proj_charge.size(); i++){
    // std::cout << proj_channel.at(i) << " " << proj_timeslice.at(i) << " " << proj_charge.at(i) << " " << proj_charge_err.at(i) << std::endl;
    
    if (proj_channel.at(i) < 2400){
      auto it = proj_data_u_map.find(std::make_pair(proj_channel.at(i),proj_timeslice.at(i)));
      if (it == proj_data_u_map.end()){
	proj_data_u_map[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_pair(proj_charge.at(i),proj_charge_err.at(i));
      }else{
	it->second.first += proj_charge.at(i);
	it->second.second = sqrt(pow(it->second.second,2) + pow(proj_charge_err.at(i),2));
      }
    }else if (proj_channel.at(i) < 4800){
      auto it = proj_data_v_map.find(std::make_pair(proj_channel.at(i),proj_timeslice.at(i)));
      if (it == proj_data_v_map.end()){
	proj_data_v_map[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_pair(proj_charge.at(i),proj_charge_err.at(i));
      }else{
	it->second.first += proj_charge.at(i);
	it->second.second = sqrt(pow(it->second.second,2) + pow(proj_charge_err.at(i),2));
      }
    }else{
      auto it = proj_data_w_map.find(std::make_pair(proj_channel.at(i),proj_timeslice.at(i)));
      if (it == proj_data_w_map.end()){
	proj_data_w_map[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_pair(proj_charge.at(i),proj_charge_err.at(i));
      }else{
	it->second.first += proj_charge.at(i);
	it->second.second = sqrt(pow(it->second.second,2) + pow(proj_charge_err.at(i),2));
      }
    }
  }

  std::cout << "dQ/dx fit: " << fine_tracking_path.size() << " " << proj_data_u_map.size() << " " << proj_data_v_map.size() << " " << proj_data_w_map.size() << std::endl;
  
  int n_3D_pos = fine_tracking_path.size();
  int n_2D_u = proj_data_u_map.size();
  int n_2D_v = proj_data_v_map.size();
  int n_2D_w = proj_data_w_map.size();

  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w), pred_data_u_2D(n_2D_u);
  
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  // initialize 
  for (int i = 0;i!=n_3D_pos;i++){
    pos_3D_init(i) = 50000;
  }

  
  // loop through the 3D points ...
  for (int i=0;i!=n_3D_pos;i++){
    

    Point prev_rec_pos, next_rec_pos;
    Point curr_rec_pos(fine_tracking_path.at(i).x, fine_tracking_path.at(i).y, fine_tracking_path.at(i).z);

    if (i==0){
      next_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i+1).x)/2.;
      next_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i+1).y)/2.;
      next_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i+1).z)/2.;
      double length = sqrt(pow(fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x,2)
			   +pow(fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y,2)
			   +pow(fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z,2) );
      prev_rec_pos.x = fine_tracking_path.at(i).x - (fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x)/length * 1.5*units::mm;
      prev_rec_pos.y = fine_tracking_path.at(i).y - (fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y)/length * 1.5*units::mm ;
      prev_rec_pos.z = fine_tracking_path.at(i).z - (fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z)/length * 1.5*units::mm;
    }else if (i==n_3D_pos-1){
      prev_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i-1).x)/2.;
      prev_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i-1).y)/2.;
      prev_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i-1).z)/2.;
      double length = sqrt(pow(fine_tracking_path.at(i-1).x-fine_tracking_path.at(i).x,2)
			   +pow(fine_tracking_path.at(i-1).y-fine_tracking_path.at(i).y,2)
			   +pow(fine_tracking_path.at(i-1).z-fine_tracking_path.at(i).z,2) );
      next_rec_pos.x = fine_tracking_path.at(i).x - (fine_tracking_path.at(i-1).x-fine_tracking_path.at(i).x);///length * 1.5*units::mm;
      next_rec_pos.y = fine_tracking_path.at(i).y - (fine_tracking_path.at(i-1).y-fine_tracking_path.at(i).y);///length * 1.5*units::mm ;
      next_rec_pos.z = fine_tracking_path.at(i).z - (fine_tracking_path.at(i-1).z-fine_tracking_path.at(i).z);///length * 1.5*units::mm;
      
    }else{
      prev_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i-1).x)/2.;
      prev_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i-1).y)/2.;
      prev_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i-1).z)/2.;

      next_rec_pos.x = (fine_tracking_path.at(i).x+fine_tracking_path.at(i+1).x)/2.;
      next_rec_pos.y = (fine_tracking_path.at(i).y+fine_tracking_path.at(i+1).y)/2.;
      next_rec_pos.z = (fine_tracking_path.at(i).z+fine_tracking_path.at(i+1).z)/2.;
    }

    // not needed ... 
    /* curr_rec_pos.x = (curr_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    /* prev_rec_pos.x = (prev_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    /* next_rec_pos.x = (next_rec_pos.x + 0.6*units::cm)/1.098 * 1.101 - 1 * 0.1101*units::cm ; */
    
    dx.push_back(sqrt(pow(curr_rec_pos.x-prev_rec_pos.x,2)+pow(curr_rec_pos.y-prev_rec_pos.y,2)+pow(curr_rec_pos.z-prev_rec_pos.z,2))
		 +sqrt(pow(curr_rec_pos.x-next_rec_pos.x,2)+pow(curr_rec_pos.y-next_rec_pos.y,2)+pow(curr_rec_pos.z-next_rec_pos.z,2)));

    
    std::vector<double> centers_U ;
    std::vector<double> centers_V ;
    std::vector<double> centers_W ;
    std::vector<double> centers_T ;
    std::vector<double> sigmas_T;
    std::vector<double> sigmas_U;
    std::vector<double> sigmas_V;
    std::vector<double> sigmas_W;
    std::vector<double> weights;


    
    for (int j=0;j!=5;j++){
      Point reco_pos;
      reco_pos.x = prev_rec_pos.x + (curr_rec_pos.x-prev_rec_pos.x)/5.*(j+0.5);
      reco_pos.y = prev_rec_pos.y + (curr_rec_pos.y-prev_rec_pos.y)/5.*(j+0.5);
      reco_pos.z = prev_rec_pos.z + (curr_rec_pos.z-prev_rec_pos.z)/5.*(j+0.5);
      double central_T = offset_t + slope_xt * reco_pos.x ;
      double central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
      double central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400;
      double central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800;
      double weight = sqrt(pow(prev_rec_pos.x-curr_rec_pos.x,2)+
			   pow(prev_rec_pos.y-curr_rec_pos.y,2)+
			   pow(prev_rec_pos.z-curr_rec_pos.z,2));

      double drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ;
      //  std::cout << drift_time/units::microsecond << std::endl;
      double diff_sigma_L = sqrt(2* DL * drift_time);
      double diff_sigma_T = sqrt(2* DT * drift_time);

      double sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width;
      double sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u;
      double sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v;
      double sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w;
      
      centers_U.push_back(central_U);
      centers_V.push_back(central_V);
      centers_W.push_back(central_W);
      centers_T.push_back(central_T);

      weights.push_back(weight);
      
      sigmas_U.push_back(sigma_T_u);
      sigmas_V.push_back(sigma_T_v);
      sigmas_W.push_back(sigma_T_w);
      sigmas_T.push_back(sigma_L);


      reco_pos.x = next_rec_pos.x + (curr_rec_pos.x-next_rec_pos.x)/5.*(j+0.5);
      reco_pos.y = next_rec_pos.y + (curr_rec_pos.y-next_rec_pos.y)/5.*(j+0.5);
      reco_pos.z = next_rec_pos.z + (curr_rec_pos.z-next_rec_pos.z)/5.*(j+0.5);
      central_T = offset_t + slope_xt * reco_pos.x ;
      central_U = offset_u + (slope_yu * reco_pos.y + slope_zu * reco_pos.z);
      central_V = offset_v + (slope_yv * reco_pos.y + slope_zv * reco_pos.z)+2400;
      central_W = offset_w + (slope_yw * reco_pos.y + slope_zw * reco_pos.z)+4800;
      weight = sqrt(pow(next_rec_pos.x-curr_rec_pos.x,2)+
		    pow(next_rec_pos.y-curr_rec_pos.y,2)+
		    pow(next_rec_pos.z-curr_rec_pos.z,2));
      
      drift_time = reco_pos.x/time_slice_width * nrebin * 0.5*units::microsecond  - flash_time ;
      //  std::cout << drift_time/units::microsecond << std::endl;
      diff_sigma_L = sqrt(2* DL * drift_time);
      diff_sigma_T = sqrt(2* DT * drift_time);

      sigma_L = sqrt(pow(diff_sigma_L,2) + pow(add_sigma_L,2))/time_slice_width;
      sigma_T_u = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_u_T,2))/pitch_u;
      sigma_T_v = sqrt(pow(diff_sigma_T,2) + pow(ind_sigma_v_T,2))/pitch_v;
      sigma_T_w = sqrt(pow(diff_sigma_T,2) + pow(col_sigma_w_T,2))/pitch_w;
      
      
      centers_U.push_back(central_U);
      centers_V.push_back(central_V);
      centers_W.push_back(central_W);
      centers_T.push_back(central_T);

      weights.push_back(weight);
      
      sigmas_U.push_back(sigma_T_u);
      sigmas_V.push_back(sigma_T_v);
      sigmas_W.push_back(sigma_T_w);
      sigmas_T.push_back(sigma_L);
    }

    // fill in the data ... 
    if (i==0){
      int n_u = 0;
      for (auto it = proj_data_u_map.begin(); it!= proj_data_u_map.end(); it++){
	data_u_2D(n_u) = it->second.first/it->second.second;
	n_u ++;
      }
      int n_v = 0;
      for (auto it = proj_data_v_map.begin(); it!= proj_data_v_map.end(); it++){
	data_v_2D(n_v) = it->second.first/it->second.second;
	n_v ++;
      }
      int n_w = 0;
      for (auto it = proj_data_w_map.begin(); it!= proj_data_w_map.end(); it++){
	data_w_2D(n_w) = it->second.first/it->second.second;
	n_w ++;
      }
    }

    int n_u = 0;
    for (auto it = proj_data_u_map.begin(); it!= proj_data_u_map.end(); it++){
      if (fabs(it->first.first - centers_U.front()) <= 10 &&
      	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_U, sigmas_U, weights , 0 , 4);
	if (value > 0)
	  RU.insert(n_u,i) = value/it->second.second;
      }
      n_u ++;
    }
    
    int n_v = 0;
    for (auto it = proj_data_v_map.begin(); it!= proj_data_v_map.end(); it++){
      if (fabs(it->first.first - centers_V.front()) <= 10 &&
	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_V, sigmas_V, weights , 0 , 4);
	if (value > 0)
	  RV.insert(n_v,i) = value/it->second.second;
      }
      n_v ++;
    }
    
    int n_w = 0;

    // double sum1 = 0;
    
    for (auto it = proj_data_w_map.begin(); it!= proj_data_w_map.end(); it++){
      if (fabs(it->first.first - centers_W.front()) <= 10 &&
      	  fabs(it->first.second - centers_T.front()) <= 10 ){
	double value = cal_gaus_integral_seg(it->first.second, it->first.first,centers_T, sigmas_T, centers_W, sigmas_W, weights , 0 , 4);
	//	sum1 += value;
	/* if (fabs(it->first.first - centers_W.front()) > 10 ||  */
      /* 	  fabs(it->first.second - centers_T.front()) > 10 ) */
      /* 	if (value!=0) */
      //      std::cout << value << " haha" << std::endl;  
	if (value>0)
	  RW.insert(n_w,i) = value/it->second.second;
      }
      n_w ++;
    }
    // std::cout << sum1 << std::endl;
    
  }
  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());

  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos) ;
  for (size_t i=0;i!=n_3D_pos;i++){
    if (i==0){
      FMatrix.insert(0,0) = -1;
      FMatrix.insert(0,1) = 1.;
    }else if (i==n_3D_pos-1){
      FMatrix.insert(i,i) = -1.;
      FMatrix.insert(i,i-1) = 1.;
    }else{
      FMatrix.insert(i,i)=-2.;
      FMatrix.insert(i,i+1)=1.;
      FMatrix.insert(i,i-1)=1.;
    }
  }
  double lambda = 0.0;
  FMatrix *= lambda;
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;
  Eigen::SparseMatrix<double> A =  RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;//
  solver.compute(A);
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);
  if (std::isnan(solver.error())){
    pos_3D = solver.solve(b);
  }

  double sum = 0 ;
  for (int i=0;i!=n_3D_pos;i++){
    std::cout << i << " "<< pos_3D(i) << " " << dx.at(i)/units::cm << " " << pos_3D(i)/dx.at(i)*units::cm << std::endl;
    dQ.push_back(pos_3D(i));
    sum += pos_3D(i);
  }
  std::cout << "total: " << sum << std::endl;

  /* pred_data_u_2D = RU * pos_3D; */
  /* for (int i=0;i!=n_2D_u;i++){ */
  /*   std::cout << pred_data_u_2D(i) << " " << data_u_2D(i) << std::endl; */
  /* } */
  
  /* sum = 0; */
  /* for (int i=0;i!=n_2D_u;i++){ */
  /*   sum += data_u_2D(i); */
  /* } */
  /* std::cout << sum << std::endl; */
  /* sum = 0; */
  /* for (int i=0;i!=n_2D_v;i++){ */
  /*   sum += data_v_2D(i); */
  /* } */
  /* std::cout << sum << std::endl; */
  /* sum = 0; */
  /* for (int i=0;i!=n_2D_w;i++){ */
  /*   sum += data_w_2D(i); */
  /* } */
  /* std::cout << sum << std::endl; */
  
  //
  
  // input: fine tracking trajectory ... 
  // See the collected 2D projection points ...
  // dimension of the 
  // for (size_t i=0;i!=fine_tracking_path.size(); i++){
    //std::cout << fine_tracking_path.at(i).x/units::cm << " " << fine_tracking_path.at(i).y/units::cm << " " << fine_tracking_path.at(i).z/units::cm << std::endl;
  //}

  // Loop the data
  /* for (auto it = proj_data_u_map.begin(); it!= proj_data_u_map.end(); it++){ */
  /*   std::cout << it->first.first << " " << it->first.second << " " << it->second << std::endl; */
  /* } */
    
  // Need to take into account the software filters regarding the induction and collection ...

  
}
