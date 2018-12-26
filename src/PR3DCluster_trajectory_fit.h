void PR3DCluster::organize_wcps_path(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec, std::vector<double>& distances, double low_dis_limit){
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      path_wcps_vec.push_back(*it);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	path_wcps_vec.push_back(*it);
	distances.push_back(dis);
      }
    }
  }
}

void PR3DCluster::fill_2d_charge(std::map<std::pair<int,int>,double>& map_2D_ut_charge,std::map<std::pair<int,int>,double>& map_2D_ut_charge_err, std::map<std::pair<int,int>,double>& map_2D_vt_charge,std::map<std::pair<int,int>,double>& map_2D_vt_charge_err, std::map<std::pair<int,int>,double>& map_2D_wt_charge,std::map<std::pair<int,int>,double>& map_2D_wt_charge_err){

  for (auto it=mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    WireChargeMap& wire_charge_err_map = mcell->get_wirechargeerr_map();
    for (auto it = wire_charge_map.begin(); it!= wire_charge_map.end(); it++){
      const GeomWire* wire = it->first;
      double charge = it->second;
      double charge_err = wire_charge_err_map[wire];
      //std::cout << wire << " " << charge << " " << charge_err << std::endl;

      // hack the charge ... 
      if (charge <=0){
	continue;
	//charge = 1000;
	//charge_err = 1000;
      }
      
      if (wire->iplane()==0){
	map_2D_ut_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_ut_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }else if (wire->iplane()==1){
	map_2D_vt_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_vt_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }else{
	map_2D_wt_charge[std::make_pair(wire->index(),time_slice)] = charge;
	map_2D_wt_charge_err[std::make_pair(wire->index(),time_slice)] = charge_err;
      }
	
    }
    //    std::cout << wire_charge_map.size() << " " << wire_charge_err_map.size() << std::endl;
  }
}

void PR3DCluster::form_map_graph_based(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec, std::vector<double>& distances, std::map<int,std::set<std::pair<int,int>>>& map_3D_2DU_set, std::map<int,std::set<std::pair<int,int>>>& map_3D_2DV_set, std::map<int,std::set<std::pair<int,int>>>& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set){

  WireCell::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  
  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();
  
  float coef1 = 2 * pow(sin(angle_u),2);
  float coef2 = 2 * (pow(sin(angle_u),2) - pow(cos(angle_u),2));
  
  typedef boost::property_map<MCUGraph, boost::vertex_index_t>::type IndexMap;
  IndexMap index = get(boost::vertex_index,*graph);
  typedef boost::graph_traits<MCUGraph>::adjacency_iterator adjacency_iterator;
  int nlevel = 3;
  double dis_cut,time_cut;
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    int current_index = path_wcps_vec.at(i).index;
    if (i==0){
      dis_cut = std::min(distances.at(i) * 0.6,0.8*units::cm);
    }else if (i==path_wcps_vec.size()-1){
      dis_cut = std::min(distances.back() * 0.6,0.8*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*0.9,distances.at(i)*0.9),1.2*units::cm);
    }
    
    time_cut = 3; // allow +- 3 time slices and then distance cut ... 
    std::set<std::pair<int,int>> T2DU_set;
    std::set<std::pair<int,int>> T2DV_set;
    std::set<std::pair<int,int>> T2DW_set;
    map_3D_2DU_set[i] = T2DU_set;
    map_3D_2DV_set[i] = T2DV_set;
    map_3D_2DW_set[i] = T2DW_set;
    
    
    
    std::set<int> total_vertices_found;
    std::set<int> vertices_to_be_examined;
    std::set<int> vertices_saved_for_next;
    total_vertices_found.insert(current_index);
    vertices_to_be_examined.insert(current_index);

    for (int j=0;j!=nlevel;j++){
      for (auto it = vertices_to_be_examined.begin(); it!=vertices_to_be_examined.end(); it++){
     	int temp_current_index = (*it);
     	std::pair<adjacency_iterator, adjacency_iterator> neighbors = boost::adjacent_vertices(vertex(temp_current_index,*graph),*graph);
     	for (; neighbors.first!=neighbors.second; ++neighbors.first){
    	  //std::cout << *neighbors.first << " " << *neighbors.second << std::endl;
    	  if (total_vertices_found.find(index(*neighbors.first))==total_vertices_found.end()){
	    total_vertices_found.insert(index(*neighbors.first));
	    vertices_saved_for_next.insert(index(*neighbors.first));
    	  }
    	}
      }
      vertices_to_be_examined = vertices_saved_for_next;
    }
    SMGCSet nearby_mcells_set;
    for (auto it = total_vertices_found.begin(); it!=total_vertices_found.end(); it++){
      SlimMergeGeomCell *mcell = cloud.pts[*it].mcell;
      nearby_mcells_set.insert(mcell);
    }
    // std::cout << i << " " << total_vertices_found.size() << " " << nearby_mcells_set.size() << std::endl;
    

    int cur_time_slice = cloud.pts[current_index].mcell->GetTimeSlice();
    int cur_wire_u = cloud.pts[current_index].index_u;
    int cur_wire_v = cloud.pts[current_index].index_v;
    int cur_wire_w = cloud.pts[current_index].index_w;

    // if (abs(cur_time_slice-1261)==0)
    //   std::cout << "Center: " << i << " " << path_wcps_vec.at(i).mcell << " " << current_index <<  " " << cloud.pts[current_index].index << " " <<
    // 	cloud.pts[current_index].mcell << " " << cloud.pts[current_index].x << " " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_v << std::endl;

    // Now fill the other maps ...
    for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int this_time_slice = mcell->GetTimeSlice();
      double rem_dis_cut = pow(dis_cut,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      if (rem_dis_cut >0 && fabs(cur_time_slice-this_time_slice)<=time_cut){
    	//	std::cout << this_time_slice << std::endl;
    	GeomWireSelection uwires = mcell->get_uwires();
    	GeomWireSelection vwires = mcell->get_vwires();
    	GeomWireSelection wwires = mcell->get_wwires();
	
    	float min_u_dis;
    	if (cur_wire_u < uwires.front()->index()){
    	  min_u_dis = uwires.front()->index()-cur_wire_u;
    	}else if (cur_wire_u >= uwires.front()->index() &&
    		  cur_wire_u <= uwires.back()->index()){
    	  min_u_dis = 0;
    	}else{
    	  min_u_dis = cur_wire_u-uwires.back()->index();
    	}
    	float min_v_dis;
    	if (cur_wire_v < vwires.front()->index()){
    	  min_v_dis = vwires.front()->index()-cur_wire_v;
    	}else if (cur_wire_v >= vwires.front()->index() &&
    		  cur_wire_v <= vwires.back()->index()){
    	  min_v_dis = 0;
    	}else{
    	  min_v_dis = cur_wire_v-vwires.back()->index();
    	}
    	float min_w_dis;
    	if (cur_wire_w < wwires.front()->index()){
    	  min_w_dis = wwires.front()->index()-cur_wire_w;
    	}else if (cur_wire_w >= wwires.front()->index() &&
    		  cur_wire_w <= wwires.back()->index()){
    	  min_w_dis = 0;
    	}else{
    	  min_w_dis = cur_wire_w-wwires.back()->index();
    	}

	// Note this assumes w is vertical wires, U and V have equal angle between then ...

    	float range_u = rem_dis_cut*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_v = rem_dis_cut*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_w = (rem_dis_cut*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
	//std::cout << coef1 << " " << coef2 << std::endl;

    	// if (abs(cur_time_slice-1261)==0)
    	//   std::cout << min_u_dis << " " << min_v_dis << " " << min_w_dis << std::endl;
	
    	if ( range_u > 0 && range_v >0 && range_w > 0){
    	  float low_u_limit = cur_wire_u - sqrt(range_u)/pitch_u;
    	  float high_u_limit = cur_wire_u + sqrt(range_u)/pitch_u;
    	  float low_v_limit = cur_wire_v - sqrt(range_v)/pitch_v;
    	  float high_v_limit = cur_wire_v + sqrt(range_v)/pitch_v;
    	  float low_w_limit = cur_wire_w - sqrt(range_w)/pitch_w;
    	  float high_w_limit = cur_wire_w + sqrt(range_w)/pitch_w;

    	  // if (abs(cur_time_slice-1261)==0)
    	  //   std::cout << low_u_limit << " " << high_u_limit << " " << low_v_limit << " " << high_v_limit << " " << low_w_limit << " " << high_w_limit << std::endl;
	  
    	  WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    	  for (auto it1 = wire_charge_map.begin(); it1!= wire_charge_map.end(); it1++){
	    const GeomWire *wire = it1->first;
	    if (it1->second >0){
	      //if (1>0){
	      if (wire->iplane()==0){
		// U plane ...
		if (wire->index() >= low_u_limit && wire->index() <= high_u_limit){
		  map_3D_2DU_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DU_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DU_3D_set.end()){
		    std::set<int>  temp_set;
		    temp_set.insert(i);
		    map_2DU_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
		    map_2DU_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }else if (wire->iplane()==1){
		// V plane ...
		if (wire->index() >= low_v_limit && wire->index() <= high_v_limit){
		  map_3D_2DV_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DV_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DV_3D_set.end()){
		    std::set<int>  temp_set;
		    temp_set.insert(i);
		    map_2DV_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
    		  map_2DV_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }else{
		// W plane ... 
		if (wire->index() >= low_w_limit && wire->index() <= high_w_limit){
		  map_3D_2DW_set[i].insert(std::make_pair(wire->index(),this_time_slice));
		  if (map_2DW_3D_set.find(std::make_pair(wire->index(),this_time_slice))==map_2DW_3D_set.end()){
		    std::set<int> temp_set;
		    temp_set.insert(i);
		    map_2DW_3D_set[std::make_pair(wire->index(),this_time_slice)] = temp_set;
		  }else{
		    map_2DW_3D_set[std::make_pair(wire->index(),this_time_slice)].insert(i);
		  }
		}
	      }
	    }
    	  }
	}
      }
    }
  } // i loop ... 
  
  
}


void PR3DCluster::trajectory_fit(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec, std::vector<double>& distances,std::map<int,std::set<std::pair<int,int>>>& map_3D_2DU_set, std::map<int,std::set<std::pair<int,int>>>& map_3D_2DV_set, std::map<int,std::set<std::pair<int,int>>>& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>,double>& map_2D_ut_charge,std::map<std::pair<int,int>,double>& map_2D_ut_charge_err, std::map<std::pair<int,int>,double>& map_2D_vt_charge,std::map<std::pair<int,int>,double>& map_2D_vt_charge_err, std::map<std::pair<int,int>,double>& map_2D_wt_charge,std::map<std::pair<int,int>,double>& map_2D_wt_charge_err){
  
  // map 2D points to its index
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DU_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DV_index;
  std::vector<std::pair<std::pair<int,int>,int>> vec_2DW_index;

  for (auto it = map_2DU_3D_set.begin(); it!= map_2DU_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DU_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DV_3D_set.begin(); it!= map_2DV_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DV_index.push_back(std::make_pair(it->first,*it1));
    }
  }
  for (auto it = map_2DW_3D_set.begin(); it!= map_2DW_3D_set.end(); it++){
    for (auto it1 =it->second.begin(); it1!=it->second.end(); it1++){
      vec_2DW_index.push_back(std::make_pair(it->first,*it1));
    }
  }

  TPCParams& mp = Singleton<TPCParams>::Instance();
  double pitch_u = mp.get_pitch_u();
  double pitch_v = mp.get_pitch_v();
  double pitch_w = mp.get_pitch_w();
  double angle_u = mp.get_angle_u();
  double angle_v = mp.get_angle_v();
  double angle_w = mp.get_angle_w();
  double time_slice_width = mp.get_ts_width();
  double first_u_dis = mp.get_first_u_dis();
  double first_v_dis = mp.get_first_v_dis();
  double first_w_dis = mp.get_first_w_dis();


  double slope_x = 1./time_slice_width;
  double first_t_dis = path_wcps_vec.at(0).mcell->GetTimeSlice()*time_slice_width - path_wcps_vec.at(0).x;
  double offset_t = first_t_dis / time_slice_width;
  
  
  
  //  convert Z to W ... 
  double slope_zw = 1./pitch_w * cos(angle_w);
  double slope_yw = 1./pitch_w * sin(angle_w);
  
  double slope_yu = -1./pitch_u * sin(angle_u);
  double slope_zu = 1./pitch_u * cos(angle_u);
  double slope_yv = -1./pitch_v * sin(angle_v);
  double slope_zv = 1./pitch_v * cos(angle_v);
  
  //convert Y,Z to U,V
  double offset_w = -first_w_dis/pitch_w;
  double offset_u = -first_u_dis/pitch_u;
  double offset_v = -first_v_dis/pitch_v;

   // form matrix ...
  int n_3D_pos = 3 * path_wcps_vec.size();
  int n_2D_u = 2 * vec_2DU_index.size();
  int n_2D_v = 2 * vec_2DV_index.size();
  int n_2D_w = 2 * vec_2DW_index.size();
  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    pos_3D_init(3*i) = path_wcps_vec.at(i).x;
    pos_3D_init(3*i+1) = path_wcps_vec.at(i).y;
    pos_3D_init(3*i+2) = path_wcps_vec.at(i).z;
  }
  
  // fill in the measurement ...
  for (size_t index = 0; index!=vec_2DU_index.size(); index++){
    double charge = map_2D_ut_charge[vec_2DU_index.at(index).first];
    double charge_err = map_2D_ut_charge_err[vec_2DU_index.at(index).first];
    int n_divide = map_2DU_3D_set[vec_2DU_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_u_2D(2*index) =  scaling * (vec_2DU_index.at(index).first.first - offset_u);
    data_u_2D(2*index+1) = scaling * (vec_2DU_index.at(index).first.second - offset_t);

    int index_3D = vec_2DU_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RU.insert(2*index,3*index_3D+1) = scaling * slope_yu; // Y--> U
    RU.insert(2*index,3*index_3D+2) = scaling * slope_zu; // Z--> U
    RU.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T

    // if (std::isnan(scaling * (vec_2DU_index.at(index).first.first - offset_u)) ||
    // 	std::isnan(scaling * (vec_2DU_index.at(index).first.second - offset_t)) ||
    // 	std::isnan(scaling * slope_yu) ||
    // 	std::isnan(scaling * slope_zu) ||
    // 	std::isnan(scaling * slope_x))
    //   std::cout << "Wrong U" << " " << charge << " " << charge_err << " " << n_divide<< std::endl;

    //std::cout << index << " " << index_3D << std::endl;
  }
  for (size_t index = 0; index!=vec_2DV_index.size(); index++){
    double charge = map_2D_vt_charge[vec_2DV_index.at(index).first];
    double charge_err = map_2D_vt_charge_err[vec_2DV_index.at(index).first];
    int n_divide = map_2DV_3D_set[vec_2DV_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_v_2D(2*index) =  scaling * (vec_2DV_index.at(index).first.first - offset_v);
    data_v_2D(2*index+1) = scaling * (vec_2DV_index.at(index).first.second - offset_t);

    int index_3D = vec_2DV_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RV.insert(2*index,3*index_3D+1) = scaling * slope_yv; // Y--> V
    RV.insert(2*index,3*index_3D+2) = scaling * slope_zv; // Z--> V
    RV.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T
    //std::cout << index << " " << index_3D << std::endl;

    // if (std::isnan(scaling * (vec_2DV_index.at(index).first.first - offset_v)) ||
    // 	std::isnan(scaling * (vec_2DV_index.at(index).first.second - offset_t)) ||
    // 	std::isnan(scaling * slope_yv) ||
    // 	std::isnan(scaling * slope_zv) ||
    // 	std::isnan(scaling * slope_x))
    //   std::cout << "Wrong V" << " " << charge << " " << charge_err << " " << n_divide << std::endl;
  }
  for (size_t index = 0; index!=vec_2DW_index.size(); index++){
    double charge = map_2D_wt_charge[vec_2DW_index.at(index).first];
    double charge_err = map_2D_wt_charge_err[vec_2DW_index.at(index).first];
    int n_divide = map_2DW_3D_set[vec_2DW_index.at(index).first].size();
    double scaling = charge/charge_err/n_divide;
    data_w_2D(2*index) =  scaling * (vec_2DW_index.at(index).first.first - offset_w);
    data_w_2D(2*index+1) = scaling * (vec_2DW_index.at(index).first.second - offset_t);

    int index_3D = vec_2DW_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    RW.insert(2*index,3*index_3D+2) = scaling * slope_zw; // Z--> W
    RW.insert(2*index+1,3*index_3D) = scaling * slope_x; // X --> T

     // if (std::isnan(scaling * (vec_2DW_index.at(index).first.first - offset_w)) ||
     // 	std::isnan(scaling * (vec_2DW_index.at(index).first.second - offset_t)) ||
     // 	 std::isnan(scaling * slope_zw) ||
     // 	std::isnan(scaling * slope_x))
     //  std::cout << "Wrong W" << " " << charge << " " << charge_err << " " << n_divide << std::endl;
  }
  
  Eigen::SparseMatrix<double> RUT = Eigen::SparseMatrix<double>(RU.transpose());
  Eigen::SparseMatrix<double> RVT = Eigen::SparseMatrix<double>(RV.transpose());
  Eigen::SparseMatrix<double> RWT = Eigen::SparseMatrix<double>(RW.transpose());

   double lambda = 2 // strength 
    * sqrt(9. //  average chi2 guessted
 	   * (vec_2DU_index.size() + vec_2DV_index.size() + vec_2DW_index.size()) // how many of them
 	   * 6 * 6 // charge/charge_err estimation ... 
 	   /(path_wcps_vec.size() * 1.)); //weighting
  lambda = 0;
  //double dis_range = 0.5*units::cm/sqrt(3.);
  double angle_range = 0.25;
 
  //std::cout << lambda/ dis_range << std::endl;
  
  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos) ;
  //Eigen::SparseMatrix<double> PMatrix(n_3D_pos, n_3D_pos) ;
  // distances[i]
  // 2nd order ...
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    // PMatrix.insert(3*i,3*i)=1;
    // PMatrix.insert(3*i+1,3*i+1)=1;
    // PMatrix.insert(3*i+2,3*i+2)=1;
    if (i==0){
      FMatrix.insert(0,0) = -1./distances.at(0); // X
      FMatrix.insert(0,3) = 1./distances.at(0);
      FMatrix.insert(1,1) = -1./distances.at(0); // Y
      FMatrix.insert(1,4) = 1./distances.at(0);
      FMatrix.insert(2,2) = -1./distances.at(0); // Z
      FMatrix.insert(2,5) = 1./distances.at(0);
    }else if (i==path_wcps_vec.size()-1){
      FMatrix.insert(3*i,3*i) = -1./distances.at(path_wcps_vec.size()-2); // X
      FMatrix.insert(3*i,3*i-3) = 1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(path_wcps_vec.size()-2);
      FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(path_wcps_vec.size()-2);
    }else{
      FMatrix.insert(3*i,3*i-3) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i,3*i) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i,3*i+3) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i+1,3*i+4) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(i-1)/(distances.at(i)+distances.at(i-1));
      FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(i-1)/distances.at(i);
      FMatrix.insert(3*i+2,3*i+5) = 1./distances.at(i)/(distances.at(i)+distances.at(i-1));
    }
  }

  
  FMatrix *= lambda/angle_range ; // disable the angue cut ... 
  
  
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  // Eigen::SparseMatrix<double> PMatrixT = Eigen::SparseMatrix<double>(PMatrix.transpose());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;// + PMatrixT * pos_3D_init * pow(lambda/dis_range,2);

   Eigen::SparseMatrix<double> A =   RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;// + PMatrixT * PMatrix * pow(lambda/dis_range,2);

   solver.compute(A);
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);

   if (std::isnan(solver.error())){
    pos_3D = solver.solve(b);
    //std::cout << "#iterations: " << solver.iterations() << std::endl;
    //std::cout << "#estimated error: " << solver.error() << std::endl;
  }
  
  //std::cout << path_wcps_vec.size() << " " << map_2DU_index.size() << " " << map_2DV_index.size() << " " << map_2DW_index.size() << std::endl;

  if (std::isnan(solver.error())){
    fine_tracking_path.clear();
    if (fine_tracking_path.size()==0){
      for (size_t i=0;i!=path_wcps_vec.size();i++){
	Point p;
	p.x = path_wcps_vec.at(i).x;
	p.y = path_wcps_vec.at(i).y;
	p.z = path_wcps_vec.at(i).z;
	fine_tracking_path.push_back(p);
      }
    }
  }else{
    flag_fine_tracking = true;
    fine_tracking_path.clear();
    for (size_t i=0;i!=path_wcps_vec.size();i++){
      Point p;
      p.x = pos_3D(3*i);
      p.y = pos_3D(3*i+1);
      p.z = pos_3D(3*i+2);
      if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)){
      }else{
	fine_tracking_path.push_back(p);
      }
    }
    // std::cout << p.x << " " << p.y << " " << p.z << " " << path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).y << " " << path_wcps_vec.at(i).z << std::endl;
  }
  
}


void PR3DCluster::fine_tracking(int num_pts_cut){

  // first round tracking with graph-based solution
  if (path_wcps.size() < num_pts_cut) return;

  // organize trajectory point based on graph-based solution 
  double low_dis_limit = 0.5*units::cm;
  std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec;
  std::vector<double> distances;
  organize_wcps_path(path_wcps_vec,distances,low_dis_limit);

  // fill initial charge ... 
  std::map<std::pair<int,int>,double> map_2D_ut_charge, map_2D_ut_charge_err;
  std::map<std::pair<int,int>,double> map_2D_vt_charge, map_2D_vt_charge_err;
  std::map<std::pair<int,int>,double> map_2D_wt_charge, map_2D_wt_charge_err;
  fill_2d_charge(map_2D_ut_charge, map_2D_ut_charge_err, map_2D_vt_charge, map_2D_vt_charge_err, map_2D_wt_charge, map_2D_wt_charge_err);
  
  
  // map index ... 
  // map 3D index to set of 2D points
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DU_set;
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DV_set;
  std::map<int,std::set<std::pair<int,int>>> map_3D_2DW_set;
  // map 2D points to 3D indices
  std::map<std::pair<int,int>,std::set<int>> map_2DU_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DV_3D_set;
  std::map<std::pair<int,int>,std::set<int>> map_2DW_3D_set;
  form_map_graph_based(path_wcps_vec, distances,
		       map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
		       map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set);

  // trajectory fitting 
  trajectory_fit(path_wcps_vec, distances,
		 map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
		 map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
		 map_2D_ut_charge, map_2D_ut_charge_err, map_2D_vt_charge,
		 map_2D_vt_charge_err, map_2D_wt_charge, map_2D_wt_charge_err);
  
  
}
