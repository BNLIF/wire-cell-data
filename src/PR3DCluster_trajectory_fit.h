#include "TVector3.h"

void PR3DCluster::organize_wcps_vec_path(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec,  double low_dis_limit){
  std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec = path_wcps_vec;

  path_wcps_vec.clear();
  for (size_t i=0;i!=temp_wcps_vec.size();i++){
    if (path_wcps_vec.size()==0){
      path_wcps_vec.push_back(temp_wcps_vec.at(i));
    }else if (i+1==temp_wcps_vec.size()){
      double dis = sqrt(pow(temp_wcps_vec.at(i).x - path_wcps_vec.back().x,2)
  			+pow(temp_wcps_vec.at(i).y - path_wcps_vec.back().y,2)
  			+pow(temp_wcps_vec.at(i).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit * 0.75){
  	path_wcps_vec.push_back(temp_wcps_vec.at(i));
      }
    }else{
      double dis = sqrt(pow(temp_wcps_vec.at(i).x - path_wcps_vec.back().x,2)
  			+pow(temp_wcps_vec.at(i).y - path_wcps_vec.back().y,2)
  			+pow(temp_wcps_vec.at(i).z - path_wcps_vec.back().z,2));
      double dis1 = sqrt(pow(temp_wcps_vec.at(i+1).x - path_wcps_vec.back().x,2)
  			 +pow(temp_wcps_vec.at(i+1).y - path_wcps_vec.back().y,2)
  			 +pow(temp_wcps_vec.at(i+1).z - path_wcps_vec.back().z,2));
      
      if (dis > low_dis_limit ||
  	  (dis1 > low_dis_limit * 1.7
  	   && dis > low_dis_limit * 0.75)){
  	path_wcps_vec.push_back(temp_wcps_vec.at(i));
      }
    }
  }
  
  temp_wcps_vec = path_wcps_vec;
  std::vector<double> temp_angle_vec;
  for (size_t i=0;i!=temp_wcps_vec.size();i++){
    if (i==0 || i+1 == temp_wcps_vec.size()){
      temp_angle_vec.push_back(3.1415926);
    }else{
      TVector3 v1(temp_wcps_vec.at(i).x - temp_wcps_vec.at(i-1).x, temp_wcps_vec.at(i).y - temp_wcps_vec.at(i-1).y, temp_wcps_vec.at(i).z - temp_wcps_vec.at(i-1).z);
      TVector3 v2(temp_wcps_vec.at(i).x - temp_wcps_vec.at(i+1).x, temp_wcps_vec.at(i).y - temp_wcps_vec.at(i+1).y, temp_wcps_vec.at(i).z - temp_wcps_vec.at(i+1).z);
      temp_angle_vec.push_back(v1.Angle(v2));
    }
    /* if (cluster_id==9) */
    /*   std::cout << i << " " << temp_angle_vec.back()/3.1415926*180. << std::endl; */
  }

  path_wcps_vec.clear();
  for (size_t i=0;i!=temp_wcps_vec.size();i++){
    if (temp_angle_vec.at(i)/3.1415926*180.>= 75)
      path_wcps_vec.push_back(temp_wcps_vec.at(i));
  }

  
 

}


void PR3DCluster::organize_wcps_path(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec, double low_dis_limit){
  
  // copy list into a vector
  std::vector<WCPointCloud<double>::WCPoint> temp_wcps_vec(std::begin(path_wcps), std::end(path_wcps));

  for (size_t i=0;i!=temp_wcps_vec.size(); i++){
    //if (temp_wcps_vec.at(i).mcell==0) continue;
    
    if (path_wcps_vec.size()==0){
      path_wcps_vec.push_back(temp_wcps_vec.at(i));
    }else if (i+1==temp_wcps_vec.size()){
      double dis = sqrt(pow(temp_wcps_vec.at(i).x - path_wcps_vec.back().x,2)
			+pow(temp_wcps_vec.at(i).y - path_wcps_vec.back().y,2)
			+pow(temp_wcps_vec.at(i).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit * 0.5){
	path_wcps_vec.push_back(temp_wcps_vec.at(i));
      }
    }else{
      double dis = sqrt(pow(temp_wcps_vec.at(i).x - path_wcps_vec.back().x,2)
			+pow(temp_wcps_vec.at(i).y - path_wcps_vec.back().y,2)
			+pow(temp_wcps_vec.at(i).z - path_wcps_vec.back().z,2));
      double dis1 = sqrt(pow(temp_wcps_vec.at(i+1).x - path_wcps_vec.back().x,2)
			+pow(temp_wcps_vec.at(i+1).y - path_wcps_vec.back().y,2)
			+pow(temp_wcps_vec.at(i+1).z - path_wcps_vec.back().z,2));

      //  std::cout << dis/units::cm << std::endl;
      
      if (dis > low_dis_limit ||
	  (dis1 > low_dis_limit * 1.7
	  && dis > low_dis_limit * 0.5)){
	path_wcps_vec.push_back(temp_wcps_vec.at(i));
      }
    }
  }

  
  /* for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){ */
  /*   if (path_wcps_vec.size()==0){ */
  /*     path_wcps_vec.push_back(*it); */
  /*   }else{ */
  /*     double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2) */
  /* 			+pow((*it).y - path_wcps_vec.back().y,2) */
  /* 			+pow((*it).z - path_wcps_vec.back().z,2)); */
      
  /*     if (dis > low_dis_limit){ */
  /* 	path_wcps_vec.push_back(*it); */
  /* 	distances.push_back(dis); */
  /*     } */
  /*   } */
  /* } */

  
}

void PR3DCluster::fill_2d_charge(std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>, std::tuple<double,double, int> >& map_2D_wt_charge){

  std::set<SlimMergeGeomCell*> cluster_mcells_set;
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = *it;
    cluster_mcells_set.insert(mcell);
  }
  
  
  for (auto it=mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    int time_slice = mcell->GetTimeSlice();
    std::map<const GeomWire*, SMGCSelection >& timeslice_wc_map = global_wc_map[time_slice];

    WireChargeMap& wire_charge_map = mcell->get_wirecharge_map();
    WireChargeMap& wire_charge_err_map = mcell->get_wirechargeerr_map();
    for (auto it = wire_charge_map.begin(); it!= wire_charge_map.end(); it++){
      const GeomWire* wire = it->first;
      double charge = it->second;
      double charge_err = wire_charge_err_map[wire];
      // 

      // hack the charge ... 
      if (charge <=0){
	continue;
	//charge = 1000;
	//charge_err = 1000;
      }

      if (timeslice_wc_map[wire].size()>1){
	for (auto it1 = timeslice_wc_map[wire].begin(); it1!=timeslice_wc_map[wire].end(); it1++){
	  SlimMergeGeomCell *mcell1 = *it1;
	  if (cluster_mcells_set.find(mcell1)==cluster_mcells_set.end())
	    charge_err = 8000;
	    }
      }

      //  std::cout << wire << " " << charge << " " << charge_err << std::endl;
      
      if (wire->iplane()==0){
	map_2D_ut_charge[std::make_pair(wire->index(),time_slice)] = std::make_tuple(charge,charge_err,1);
      }else if (wire->iplane()==1){
	map_2D_vt_charge[std::make_pair(wire->index(),time_slice)] = std::make_tuple(charge,charge_err,1);
      }else{
	map_2D_wt_charge[std::make_pair(wire->index(),time_slice)] = std::make_tuple(charge,charge_err,1);
      }
	
    }
    //    std::cout << wire_charge_map.size() << " " << wire_charge_err_map.size() << std::endl;
  }
}

void PR3DCluster::form_map_projection_based(PointVector& ps_vec, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int, std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_wt_charge, double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){

  std::vector<double> distances;
  for (size_t i=0;i+1!=ps_vec.size();i++){
    distances.push_back(sqrt(pow(ps_vec.at(i+1).x-ps_vec.at(i).x,2) +
			     pow(ps_vec.at(i+1).y-ps_vec.at(i).y,2) +
			     pow(ps_vec.at(i+1).z-ps_vec.at(i).z,2)));
  }

  
  map_3D_2DU_set.clear();
  map_3D_2DV_set.clear();
  map_3D_2DW_set.clear();

  map_2DU_3D_set.clear();
  map_2DV_3D_set.clear();
  map_2DW_3D_set.clear();

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


  for (size_t i=0; i!=ps_vec.size(); i++){
    double dis_cut;
    
    if (i==0){
      dis_cut = std::min(distances.at(i) * end_point_factor,4/3.*end_point_factor*units::cm);
    }else if (i+1==ps_vec.size()){
      dis_cut = std::min(distances.back() * end_point_factor,4/3.*end_point_factor*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*mid_point_factor,distances.at(i)*mid_point_factor),4/3.*mid_point_factor*units::cm);
    }
    //    std::cout << i <<  " " << dis_cut/units::cm << std::endl;
    

    std::set<std::pair<int,int>> T2DU_set;
    std::set<std::pair<int,int>> T2DV_set;
    std::set<std::pair<int,int>> T2DW_set;
    map_3D_2DU_set[i] = std::make_pair(T2DU_set,0);
    map_3D_2DV_set[i] = std::make_pair(T2DV_set,0);
    map_3D_2DW_set[i] = std::make_pair(T2DW_set,0);
    
    double central_T = offset_t + slope_xt * ps_vec.at(i).x ; 
    double central_U = offset_u + (slope_yu * ps_vec.at(i).y + slope_zu * ps_vec.at(i).z);
    double central_V = offset_v + (slope_yv * ps_vec.at(i).y + slope_zv * ps_vec.at(i).z);
    double central_W = offset_w + (slope_yw * ps_vec.at(i).y + slope_zw * ps_vec.at(i).z);
    //    std::cout << central_T << " " << central_U << " " << central_V << " " << central_W << std::endl;

    std::set<int> temp_types_u;
    std::set<int> temp_types_v;
    std::set<int> temp_types_w;

    for (int j=central_T - time_cut-1; j<=central_T+time_cut+1;j++){
      // U plane ...
      for (int k=central_U - nlevel-1;k<=central_U+nlevel+1;k++){
	if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_U)*pitch_u,2)) <= dis_cut){
	   auto it1 = map_2D_ut_charge.find(std::make_pair(k,j));
	   if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > charge_cut)
	     temp_types_u.insert(std::get<2>(it1->second));
	}
      }

      // V plane ...
      for (int k=central_V - nlevel-1;k<=central_V+nlevel+1;k++){
	if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_V)*pitch_v,2)) <= dis_cut){
	   auto it1 = map_2D_vt_charge.find(std::make_pair(k,j));
	   if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > charge_cut)
	     temp_types_v.insert(std::get<2>(it1->second));
	}
      }

      // W plane ...
      for (int k=central_W - nlevel-1;k<=central_W + nlevel + 1;k++){
	if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_W)*pitch_w,2)) <= dis_cut){
	   auto it1 = map_2D_wt_charge.find(std::make_pair(k,j));
	   if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > charge_cut)
	     temp_types_w.insert(std::get<2>(it1->second));
	}
      }
    }


    
    
    
    for (int j=central_T - time_cut-1; j<=central_T+time_cut+1;j++){

      // U plane ...
      if (temp_types_u.find(0)!=temp_types_u.end() && temp_types_u.size()==1){
      }else{
	for (int k=central_U - nlevel-1;k<=central_U+nlevel+1;k++){
	  if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_U)*pitch_u,2)) <= dis_cut){
	    auto it1 = map_2D_ut_charge.find(std::make_pair(k,j));
	    if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > charge_cut){
	      map_3D_2DU_set[i].first.insert(std::make_pair(k,j));
	      if (temp_types_u.size()==2) map_3D_2DU_set[i].second = 1;
	      if (map_2DU_3D_set.find(std::make_pair(k,j))==map_2DU_3D_set.end()){
		std::set<int>  temp_set;
		temp_set.insert(i);
		map_2DU_3D_set[std::make_pair(k,j)] = temp_set;
	      }else{
		map_2DU_3D_set[std::make_pair(k,j)].insert(i);
	      }
	    }
	  }
	} 
      }//  U plane
      
      // V plane ...
      if (temp_types_v.find(0)!=temp_types_v.end() && temp_types_v.size()==1){
      }else{
	for (int k=central_V - nlevel-1;k<=central_V+nlevel+1;k++){
	  if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_V)*pitch_v,2)) <= dis_cut){
	    auto it1 = map_2D_vt_charge.find(std::make_pair(k,j));
	    if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > charge_cut ){
	      map_3D_2DV_set[i].first.insert(std::make_pair(k,j));
	      if (temp_types_v.size()==2) map_3D_2DV_set[i].second = 1;
	      if (map_2DV_3D_set.find(std::make_pair(k,j))==map_2DV_3D_set.end()){
		std::set<int>  temp_set;
		temp_set.insert(i);
		map_2DV_3D_set[std::make_pair(k,j)] = temp_set;
	      }else{
		map_2DV_3D_set[std::make_pair(k,j)].insert(i);
	      }
	    }
	  }
	} 
      }//  V plane
      
      // W plane ...
      if (temp_types_w.find(0)!=temp_types_w.end() && temp_types_w.size()==1){
      }else{
	for (int k=central_W - nlevel-1;k<=central_W+nlevel+1;k++){
	  if (sqrt(pow((j-central_T)*time_slice_width,2)+pow((k-central_W)*pitch_w,2)) <= dis_cut){
	    auto it1 = map_2D_wt_charge.find(std::make_pair(k,j));
	    if (it1 !=map_2D_wt_charge.end() && std::get<0>(it1->second) > charge_cut ){
	      map_3D_2DW_set[i].first.insert(std::make_pair(k,j));
	      if (temp_types_w.size()==2) map_3D_2DW_set[i].second = 1;
	      if (map_2DW_3D_set.find(std::make_pair(k,j))==map_2DW_3D_set.end()){
		std::set<int>  temp_set;
		temp_set.insert(i);
		map_2DW_3D_set[std::make_pair(k,j)] = temp_set;
	      }else{
		map_2DW_3D_set[std::make_pair(k,j)].insert(i);
	      }
	    }
	  }
	}
      } //  W plane


    }
  }

  
}

void PR3DCluster::form_map_graph_based(std::vector<WCPointCloud<double>::WCPoint>& path_wcps_vec, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>,std::tuple<double,double, int> >& map_2D_wt_charge,
				       std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int>>& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set, double end_point_factor, double mid_point_factor, int nlevel, double time_cut, double charge_cut){


  /* for (auto it = mcells.begin(); it!=mcells.end(); it++){ */
  /*   SlimMergeGeomCell *mcell = *it; */
  /*   int this_time_slice = mcell->GetTimeSlice(); */
  /*   GeomWireSelection uwires = mcell->get_uwires(); */
  /*   GeomWireSelection vwires = mcell->get_vwires(); */
  /*   GeomWireSelection wwires = mcell->get_wwires(); */
  /*   if (this_time_slice == 2387 && cluster_id == 18) */
  /*     std::cout << this_time_slice << " " << uwires.front()->index() << " " << uwires.back()->index() << " " << vwires.front()->index()+2400 << " " << vwires.back()->index()+2400 << " " << wwires.front()->index()+4800 << " " << wwires.back()->index()+4800 << std::endl; */
  /* } */
  
  std::vector<double> distances;
  for (size_t i=0;i+1!=path_wcps_vec.size();i++){
    distances.push_back(sqrt(pow(path_wcps_vec.at(i+1).x-path_wcps_vec.at(i).x,2) +
			     pow(path_wcps_vec.at(i+1).y-path_wcps_vec.at(i).y,2) +
			     pow(path_wcps_vec.at(i+1).z-path_wcps_vec.at(i).z,2)));
  }
  
  map_3D_2DU_set.clear();
  map_3D_2DV_set.clear();
  map_3D_2DW_set.clear();

  map_2DU_3D_set.clear();
  map_2DV_3D_set.clear();
  map_2DW_3D_set.clear();
  
  
  WCP::WCPointCloud<double>& cloud = point_cloud->get_cloud();
  
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
  double dis_cut;
  for (size_t i=0;i!=path_wcps_vec.size();i++){
    int current_index = path_wcps_vec.at(i).index;
    if (i==0){
      dis_cut = std::min(distances.at(i) * end_point_factor,4/3.*end_point_factor*units::cm);
    }else if (i+1==path_wcps_vec.size()){
      dis_cut = std::min(distances.back() * end_point_factor,4/3.*end_point_factor*units::cm);
    }else{
      dis_cut = std::min(std::max(distances.at(i-1)*mid_point_factor,distances.at(i)*mid_point_factor),4/3.*mid_point_factor*units::cm);
    }
    
    //  time_cut = 3; // allow +- 3 time slices and then distance cut ... 
    std::set<std::pair<int,int>> T2DU_set;
    std::set<std::pair<int,int>> T2DV_set;
    std::set<std::pair<int,int>> T2DW_set;
    map_3D_2DU_set[i] = std::make_pair(T2DU_set,0);
    map_3D_2DV_set[i] = std::make_pair(T2DV_set,0);
    map_3D_2DW_set[i] = std::make_pair(T2DW_set,0);
    
    
    
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
    
    //    std::cout << i << " " << total_vertices_found.size() << " " << nearby_mcells_set.size() << std::endl;
    

    int cur_time_slice = cloud.pts[current_index].mcell->GetTimeSlice();
    int cur_wire_u = cloud.pts[current_index].index_u;
    int cur_wire_v = cloud.pts[current_index].index_v;
    int cur_wire_w = cloud.pts[current_index].index_w;

    /* if (cluster_id==18 && i==6) */
    /*   std::cout << i << " " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v+2400 << " " << cur_wire_w+4800 << std::endl; */
    
    // if (abs(cur_time_slice-1261)==0)
    //   std::cout << "Center: " << i << " " << path_wcps_vec.at(i).mcell << " " << current_index <<  " " << cloud.pts[current_index].index << " " <<
    // 	cloud.pts[current_index].mcell << " " << cloud.pts[current_index].x << " " << cur_time_slice << " " << cur_wire_u << " " << cur_wire_v << " " << cur_wire_v << std::endl;

    
    std::set<int> temp_types_u;
    std::set<int> temp_types_v;
    std::set<int> temp_types_w;


      /* if (i==138&&cluster_id==2) */
    /*   std::cout << temp_types_w.size() << " " << *temp_types_w.begin() << std::endl; */

    double dis_cut_u = dis_cut;
    double dis_cut_v = dis_cut;
    double dis_cut_w = dis_cut;
    
    /* if (i==6 && cluster_id==18){ */
    double max_time_slice_u = 0;
    double max_time_slice_v = 0;
    double max_time_slice_w = 0;
    for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int this_time_slice = mcell->GetTimeSlice();
      GeomWireSelection uwires = mcell->get_uwires();
      GeomWireSelection vwires = mcell->get_vwires();
    	GeomWireSelection wwires = mcell->get_wwires();
    	for (auto it1 = uwires.begin(); it1!=uwires.end(); it1++){
    	  if (fabs((*it1)->index()-cur_wire_u)<=1)
    	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_u)
    	      max_time_slice_u = fabs(this_time_slice-cur_time_slice);
    	}
    	for (auto it1 = vwires.begin(); it1!=vwires.end(); it1++){
    	  if (fabs((*it1)->index()-cur_wire_v)<=1)
    	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_v)
    	      max_time_slice_v = fabs(this_time_slice-cur_time_slice);
    	}
    	for (auto it1 = wwires.begin(); it1!=wwires.end(); it1++){
    	  if (fabs((*it1)->index()-cur_wire_w)<=1)
    	    if (fabs(this_time_slice-cur_time_slice)>max_time_slice_w)
    	      max_time_slice_w = fabs(this_time_slice-cur_time_slice);
    	}
	
    }
    //std::cout << max_time_slice_u << " " << max_time_slice_v << " " << max_time_slice_w << " " << dis_cut/time_slice_width << std::endl;
    if (max_time_slice_u * time_slice_width*1.2 < dis_cut_u)
      dis_cut_u = max_time_slice_u * time_slice_width*1.2;
    if (max_time_slice_v * time_slice_width*1.2 < dis_cut_v)
      dis_cut_v = max_time_slice_v * time_slice_width*1.2;
    if (max_time_slice_w * time_slice_width*1.2 < dis_cut_w)
      dis_cut_w = max_time_slice_w * time_slice_width*1.2;
    /* } */
    
    // Now fill the other maps ...
    for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int this_time_slice = mcell->GetTimeSlice();

     
      
      double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      if ((rem_dis_cut_u>0 || rem_dis_cut_v >0 || rem_dis_cut_w > 0 ) && fabs(cur_time_slice-this_time_slice)<=time_cut){
    	//	std::cout << this_time_slice << std::endl;
    	GeomWireSelection uwires = mcell->get_uwires();
    	GeomWireSelection vwires = mcell->get_vwires();
    	GeomWireSelection wwires = mcell->get_wwires();
	float min_u_dis;

	/* if (i==25 && cluster_id==18) */
	/*   std::cout << this_time_slice << " " << uwires.front()->index() << " " << uwires.back()->index() << " " << vwires.front()->index()+2400 << " " << vwires.back()->index()+2400 << " " << wwires.front()->index()+4800 << " " << wwires.back()->index()+4800 << std::endl; */
	
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

    	float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
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

	  /* if (i==138&&cluster_id==2) */
	  /*   std::cout << this_time_slice << " " << low_u_limit << " " << high_u_limit << " " << low_v_limit << " " << high_v_limit << " " << low_w_limit << " " << high_w_limit << std::endl; */
	  
	  for (int j = std::round(low_u_limit); j<= std::round(high_u_limit); j++){
	    auto it1 = map_2D_ut_charge.find(std::make_pair(j,this_time_slice));
	    if (it1!=map_2D_ut_charge.end() && std::get<0>(it1->second) > 0 ){
	      temp_types_u.insert(std::get<2>(it1->second));
	    }
	  }
	  
	  /* for (auto it1 = map_2D_ut_charge.begin(); it1 != map_2D_ut_charge.end(); it1++){ */
	  /*   if (it1->first.first >=low_u_limit && it1->first.first <= high_u_limit && this_time_slice == it1->first.second){ */
	  /*     temp_types_u.insert(std::get<2>(it1->second)); */
	  /*   } */
	  /* } */
	  for (int j=std::round(low_v_limit); j<=std::round(high_v_limit); j++){
	    auto it1 = map_2D_vt_charge.find(std::make_pair(j,this_time_slice));
	    if (it1!=map_2D_vt_charge.end() && std::get<0>(it1->second) > 0){
	      temp_types_v.insert(std::get<2>(it1->second));
	    }
	  }
	  
	  /* for (auto it1 = map_2D_vt_charge.begin(); it1 != map_2D_vt_charge.end(); it1++){ */
	  /*   if (it1->first.first >=low_v_limit && it1->first.first <= high_v_limit && this_time_slice == it1->first.second){ */
	  /*     temp_types_v.insert(std::get<2>(it1->second)); */
	  /*   } */
	  /* } */

	  for (int j=std::round(low_w_limit); j<=std::round(high_w_limit); j++){
	    auto it1 = map_2D_wt_charge.find(std::make_pair(j,this_time_slice));
	    
	    if (it1!=map_2D_wt_charge.end() && std::get<0>(it1->second) > 0){
	      /* if (i==138&&cluster_id==2) */
	      /* 	std::cout << it1->first.first << " " << it1->first.second << " " << std::get<0>(it1->second) << " " << std::get<1>(it1->second) << " " << std::get<2>(it1->second) << std::endl; */
	      temp_types_w.insert(std::get<2>(it1->second));
	    }
	  }
	  
	  /* for (auto it1 = map_2D_wt_charge.begin(); it1 != map_2D_wt_charge.end(); it1++){ */
	  /*   if (it1->first.first >=low_w_limit && it1->first.first <= high_w_limit && this_time_slice == it1->first.second){ */
	  /*     temp_types_w.insert(std::get<2>(it1->second)); */
	  /*   } */
	  /* } */
	}
      }
    }

  
    
    // Now fill the other maps ...
    for (auto it = nearby_mcells_set.begin(); it!=nearby_mcells_set.end(); it++){
      SlimMergeGeomCell *mcell = *it;
      int this_time_slice = mcell->GetTimeSlice();
      double rem_dis_cut_u = pow(dis_cut_u,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      double rem_dis_cut_v = pow(dis_cut_v,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      double rem_dis_cut_w = pow(dis_cut_w,2) - pow((cur_time_slice - this_time_slice)*time_slice_width,2);
      if ((rem_dis_cut_u >0 || rem_dis_cut_v > 0 || rem_dis_cut_w >0) && fabs(cur_time_slice-this_time_slice)<=time_cut){
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
    	float range_u = rem_dis_cut_u*coef1 - pow(min_v_dis*pitch_v,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_v = rem_dis_cut_v*coef1 - pow(min_u_dis*pitch_u,2) - coef2*pow(min_w_dis*pitch_w,2);
    	float range_w = (rem_dis_cut_w*coef1 - pow(min_u_dis*pitch_u,2) - pow(min_v_dis*pitch_v,2))/coef2;
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

	  // all dead or overlap cases
	  if (temp_types_u.find(0)!=temp_types_u.end() && temp_types_u.size()==1){
	  }else{
	    for (int j=std::round(low_u_limit); j<=std::round(high_u_limit);j++){
	      auto it1 = map_2D_ut_charge.find(std::make_pair(j,this_time_slice));
	      if (it1!=map_2D_ut_charge.end()&& std::get<0>(it1->second) > charge_cut){
		map_3D_2DU_set[i].first.insert(std::make_pair(it1->first.first,it1->first.second)); 
		if (temp_types_u.size()==2 ) map_3D_2DU_set[i].second = 1;
		if (map_2DU_3D_set.find(std::make_pair(it1->first.first, it1->first.second))==map_2DU_3D_set.end()){
		  std::set<int>  temp_set;
		  temp_set.insert(i);
		  map_2DU_3D_set[std::make_pair(it1->first.first, it1->first.second)] = temp_set;
		}else{
		  map_2DU_3D_set[std::make_pair(it1->first.first, it1->first.second)].insert(i);
		}
	      }
	    }
	  }
	 
	  if (temp_types_v.find(0)!=temp_types_v.end() && temp_types_v.size()==1){
	  }else{
	    for (int j=std::round(low_v_limit); j<=std::round(high_v_limit);j++){
	      auto it1 = map_2D_vt_charge.find(std::make_pair(j,this_time_slice));
	      if (it1!=map_2D_vt_charge.end()&& std::get<0>(it1->second) > charge_cut){
		map_3D_2DV_set[i].first.insert(std::make_pair(it1->first.first,it1->first.second));
		if (temp_types_v.size()==2) map_3D_2DV_set[i].second = 1;
		if (map_2DV_3D_set.find(std::make_pair(it1->first.first, it1->first.second))==map_2DV_3D_set.end()){
		  std::set<int>  temp_set;
		  temp_set.insert(i);
		  map_2DV_3D_set[std::make_pair(it1->first.first, it1->first.second)] = temp_set;
		}else{
		  map_2DV_3D_set[std::make_pair(it1->first.first, it1->first.second)].insert(i);
		}
	      }
	    }
	  }
	  

	  if (temp_types_w.find(0)!=temp_types_w.end() && temp_types_w.size()==1){
	  }else{
	    for (int j=std::round(low_w_limit); j<=std::round(high_w_limit);j++){
	      auto it1 = map_2D_wt_charge.find(std::make_pair(j,this_time_slice));
	      if (it1!=map_2D_wt_charge.end()&& std::get<0>(it1->second) > charge_cut){
		map_3D_2DW_set[i].first.insert(std::make_pair(it1->first.first,it1->first.second)); 
		if (temp_types_w.size()==2) map_3D_2DW_set[i].second = 1;
		if (map_2DW_3D_set.find(std::make_pair(it1->first.first, it1->first.second))==map_2DW_3D_set.end()){
		  std::set<int>  temp_set;
		  temp_set.insert(i);
		  map_2DW_3D_set[std::make_pair(it1->first.first, it1->first.second)] = temp_set;
		}else{
		  map_2DW_3D_set[std::make_pair(it1->first.first, it1->first.second)].insert(i);
		}
	      }
	    }
	  }
	}
      }
    }

    
  } // i loop ...

  /* if (cluster_id==18){ */
  /*   int index = 6; */
  /*   for (auto it = map_3D_2DU_set[index].begin() ; it!=map_3D_2DU_set[index].end();it++){ */
  /*     auto it1 = map_2D_ut_charge.find(std::make_pair(it->first,it->second)); */
  /*     std::cout << "U: " << it->first << " " << it->second << " " << std::get<0>(it1->second) << " " << std::get<1>(it1->second) << std::endl; */
  /*   } */
  /*   for (auto it = map_3D_2DV_set[index].begin() ; it!=map_3D_2DV_set[index].end();it++){ */
  /*     auto it1 = map_2D_vt_charge.find(std::make_pair(it->first,it->second)); */
  /*     std::cout << "V: " << it->first+2400 << " " << it->second << " " << std::get<0>(it1->second) << " " << std::get<1>(it1->second) << std::endl; */
  /*   } */
  /*   for (auto it = map_3D_2DW_set[index].begin() ; it!=map_3D_2DW_set[index].end();it++){ */
  /*      auto it1 = map_2D_wt_charge.find(std::make_pair(it->first,it->second)); */
  /*     std::cout << "W: " << it->first+4800 << " " << it->second << " " << std::get<0>(it1->second) << " " << std::get<1>(it1->second) << std::endl; */
  /*   } */
  /*   std::cout << std::endl << std::endl; */
  /* } */
  
}

void PR3DCluster::trajectory_fit(PointVector& ps_vec,std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set, std::map<std::pair<int,int>,std::set<int>>& map_2DU_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DV_3D_set, std::map<std::pair<int,int>,std::set<int>>& map_2DW_3D_set,std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_ut_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_vt_charge, std::map<std::pair<int,int>, std::tuple<double, double, int> >& map_2D_wt_charge){
  

  std::vector<double> distances;
  for (size_t i=0;i+1!=ps_vec.size();i++){
    distances.push_back(sqrt(pow(ps_vec.at(i+1).x-ps_vec.at(i).x,2) +
			     pow(ps_vec.at(i+1).y-ps_vec.at(i).y,2) +
			     pow(ps_vec.at(i+1).z-ps_vec.at(i).z,2)));
  }

  //  std::cout << "haha " << std::endl;
  
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
  double first_t_dis = path_wcps.front().mcell->GetTimeSlice()*time_slice_width - path_wcps.front().x;
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
  int n_3D_pos = 3 * ps_vec.size();
  int n_2D_u = 2 * vec_2DU_index.size();
  int n_2D_v = 2 * vec_2DV_index.size();
  int n_2D_w = 2 * vec_2DW_index.size();
  Eigen::VectorXd pos_3D(n_3D_pos), data_u_2D(n_2D_u), data_v_2D(n_2D_v), data_w_2D(n_2D_w);
  Eigen::SparseMatrix<double> RU(n_2D_u, n_3D_pos) ;
  Eigen::SparseMatrix<double> RV(n_2D_v, n_3D_pos) ;
  Eigen::SparseMatrix<double> RW(n_2D_w, n_3D_pos) ;
  Eigen::VectorXd pos_3D_init(n_3D_pos);
  for (size_t i=0;i!=ps_vec.size();i++){
    pos_3D_init(3*i) = ps_vec.at(i).x;
    pos_3D_init(3*i+1) = ps_vec.at(i).y;
    pos_3D_init(3*i+2) = ps_vec.at(i).z;
  }

  //  std::cout << "haha1 " << " " << vec_2DU_index.size() << std::endl;

  double scale_ind_to_col = 4.0;
  
  // fill in the measurement ...
  for (size_t index = 0; index!=vec_2DU_index.size(); index++){
    double charge = std::get<0>(map_2D_ut_charge[vec_2DU_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_ut_charge[vec_2DU_index.at(index).first]);
    int n_divide = map_2DU_3D_set[vec_2DU_index.at(index).first].size();
    // induction plane error is large, so they have less weights in the fit to decide X position
    // for trajectory fit, we can scale it back ...
    double scaling = charge/charge_err/n_divide*scale_ind_to_col;
    /* if (cluster_id==34)  */
    //    std::cout << index << " " << map_3D_2DU_set[index].second << " " << map_3D_2DV_set[index].second " " << charge << std::endl; 

    int index_3D = vec_2DU_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z    
    if (map_3D_2DU_set[index_3D].second==1)       scaling /= scale_ind_to_col;


    
    data_u_2D(2*index) =  scaling * (vec_2DU_index.at(index).first.first - offset_u);
    data_u_2D(2*index+1) = scaling * (vec_2DU_index.at(index).first.second - offset_t);

    //    std::cout << index << " " << index_3D << " " << charge << " " << data_u_2D(2*index) << " " << data_u_2D(2*index+1) << std::endl;
    
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

  //  std::cout << "haha11 " << " " << vec_2DV_index.size() << std::endl;
  
  for (size_t index = 0; index!=vec_2DV_index.size(); index++){
    double charge = std::get<0>(map_2D_vt_charge[vec_2DV_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_vt_charge[vec_2DV_index.at(index).first]);
    int n_divide = map_2DV_3D_set[vec_2DV_index.at(index).first].size();

    int index_3D = vec_2DV_index.at(index).second; // 3*index_3D -->x  3*index_3D+1 --> y 3*index_3D+2 --> z
    double scaling = charge/charge_err/n_divide*scale_ind_to_col;
    if (map_3D_2DV_set[index_3D].second == 1)       scaling /= scale_ind_to_col;
    
    data_v_2D(2*index) =  scaling * (vec_2DV_index.at(index).first.first - offset_v);
    data_v_2D(2*index+1) = scaling * (vec_2DV_index.at(index).first.second - offset_t);


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

  //  std::cout << "haha12 " << " " << vec_2DW_index.size() << std::endl;
  
  for (size_t index = 0; index!=vec_2DW_index.size(); index++){
    double charge = std::get<0>(map_2D_wt_charge[vec_2DW_index.at(index).first]);
    double charge_err = std::get<1>(map_2D_wt_charge[vec_2DW_index.at(index).first]);
    int n_divide = map_2DW_3D_set[vec_2DW_index.at(index).first].size();

    double scaling  = charge/charge_err/n_divide;
    
        
    
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


  //  std::cout << "haha2 " << n_3D_pos << " " << ps_vec.size() << std::endl;
 
  //std::cout << lambda/ dis_range << std::endl;
  
  Eigen::SparseMatrix<double> FMatrix(n_3D_pos, n_3D_pos) ;
  //Eigen::SparseMatrix<double> PMatrix(n_3D_pos, n_3D_pos) ;
  // distances[i]
  // 2nd order ...
  for (size_t i=0;i!=ps_vec.size();i++){
    
    if (i==0){
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
	FMatrix.insert(0,0) = -1./distances.at(0); // X
	FMatrix.insert(0,3) = 1./distances.at(0);
	
	FMatrix.insert(1,1) = -1./distances.at(0); // Y
	FMatrix.insert(1,4) = 1./distances.at(0);
	
	FMatrix.insert(2,2) = -1./distances.at(0); // Z
	FMatrix.insert(2,5) = 1./distances.at(0);
      }
    }else if (i+1==ps_vec.size()){
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
	FMatrix.insert(3*i,3*i) = -1./distances.at(ps_vec.size()-2); // X
	FMatrix.insert(3*i,3*i-3) = 1./distances.at(ps_vec.size()-2);

	FMatrix.insert(3*i+1,3*i+1) = -1./distances.at(ps_vec.size()-2);
	FMatrix.insert(3*i+1,3*i-2) = 1./distances.at(ps_vec.size()-2);

	FMatrix.insert(3*i+2,3*i+2) = -1./distances.at(ps_vec.size()-2);
	FMatrix.insert(3*i+2,3*i-1) = 1./distances.at(ps_vec.size()-2);
      }
    }else{
      if (map_3D_2DU_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ||
	  map_3D_2DU_set[i].first.size()==0 && map_3D_2DW_set[i].first.size()==0 ||
	  map_3D_2DW_set[i].first.size()==0 && map_3D_2DV_set[i].first.size()==0 ){
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
  }

  double lambda = 0;
  lambda = 2 // strength 
    * sqrt(9. //  average chi2 guessted
 	   * (vec_2DU_index.size() + vec_2DV_index.size() + vec_2DW_index.size()) // how many of them
 	   * 6 * 6 // charge/charge_err estimation ... 
 	   /(ps_vec.size() * 1.)); //weighting
  double angle_range = 0.25;
  FMatrix *= lambda/angle_range ; // disable the angle cut ... 
  
  
  Eigen::SparseMatrix<double> FMatrixT = Eigen::SparseMatrix<double>(FMatrix.transpose());
  // Eigen::SparseMatrix<double> PMatrixT = Eigen::SparseMatrix<double>(PMatrix.transpose());
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  Eigen::VectorXd b = RUT * data_u_2D + RVT * data_v_2D + RWT * data_w_2D;// + PMatrixT * pos_3D_init * pow(lambda/dis_range,2);

  Eigen::SparseMatrix<double> A =   RUT * RU + RVT * RV + RWT * RW + FMatrixT * FMatrix;// + PMatrixT * PMatrix * pow(lambda/dis_range,2);

  //  std::cout << "Solve1 " << std::endl;
  
  solver.compute(A);

  //  std::cout << "Solve " << std::endl;
  
  pos_3D = solver.solveWithGuess(b,pos_3D_init);
  
  if (std::isnan(solver.error())){
    pos_3D = solver.solve(b);
    //std::cout << "#iterations: " << solver.iterations() << std::endl;
    //std::cout << "#estimated error: " << solver.error() << std::endl;
  }
  
  //std::cout << path_wcps_vec.size() << " " << map_2DU_index.size() << " " << map_2DV_index.size() << " " << map_2DW_index.size() << std::endl;

  if (std::isnan(solver.error())){
    //std::cout << "Bad fit" << std::endl;
    fine_tracking_path.clear();
    if (fine_tracking_path.size()==0){
      for (size_t i=0;i!=ps_vec.size();i++){
	Point p;
	p.x = ps_vec.at(i).x;
	p.y = ps_vec.at(i).y;
	p.z = ps_vec.at(i).z;
	fine_tracking_path.push_back(p);
      }
    }
  }else{
    // std::cout << "Good fit" << std::endl;
    flag_fine_tracking = true;
    fine_tracking_path.clear();
    for (size_t i=0;i!=ps_vec.size();i++){
      Point p;
      p.x = pos_3D(3*i);
      p.y = pos_3D(3*i+1);
      p.z = pos_3D(3*i+2);
      if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)){
      }else{
	fine_tracking_path.push_back(p);
      }
    }

    /* // if (cluster_id==5){ */
    /*   for (size_t i=0;i!=fine_tracking_path.size();i++){ */
    /* 	double cur_t = slope_x * fine_tracking_path.at(i).x + offset_t; */
    /* 	double cur_u = slope_yu *fine_tracking_path.at(i).y + slope_zu * fine_tracking_path.at(i).z + offset_u; */
    /* 	double cur_v = slope_yv *fine_tracking_path.at(i).y + slope_zv * fine_tracking_path.at(i).z + offset_v; */
    /* 	double cur_w = slope_yw *fine_tracking_path.at(i).y + slope_zw * fine_tracking_path.at(i).z + offset_w; */

    /* 	// calculate average ...  */
    /* 	double ave_ut=0, ave_uc=0, ave_us=0; */
    /* 	double ave_vt=0, ave_vc=0, ave_vs=0; */
    /* 	double ave_wt=0, ave_wc=0, ave_ws=0; */

    /* 	for (auto it = map_3D_2DU_set[i].first.begin(); it!= map_3D_2DU_set[i].first.end(); it++){ */
    /* 	  double charge = std::get<0>(map_2D_ut_charge[*it]); */
    /* 	  double charge_err = std::get<1>(map_2D_ut_charge[*it]); */
    /* 	  int n_divide =  map_2DU_3D_set[*it].size(); */
    /* 	  double scaling = charge/charge_err/n_divide*scale_ind_to_col; */
    /* 	  if (map_3D_2DU_set[i].second==1)       scaling /= scale_ind_to_col; */

    /* 	  ave_ut += it->second * scaling; */
    /* 	  ave_uc += it->first * scaling; */
    /* 	  ave_us += scaling; */
    /* 	  //	  std::cout << i << " " << charge << " " << charge_err << " " << n_divide << std::endl; */
    /* 	} */

    /* 	for (auto it = map_3D_2DV_set[i].first.begin(); it!= map_3D_2DV_set[i].first.end(); it++){ */
    /* 	  double charge = std::get<0>(map_2D_vt_charge[*it]); */
    /* 	  double charge_err = std::get<1>(map_2D_vt_charge[*it]); */
    /* 	  int n_divide =  map_2DV_3D_set[*it].size(); */
    /* 	  double scaling = charge/charge_err/n_divide*scale_ind_to_col; */
    /* 	  if (map_3D_2DV_set[i].second==1)       scaling /= scale_ind_to_col; */

    /* 	  ave_vt += it->second * scaling; */
    /* 	  ave_vc += it->first * scaling; */
    /* 	  ave_vs += scaling; */
    /* 	  //	  std::cout << i << " " << charge << " " << charge_err << " " << n_divide << std::endl; */
    /* 	} */

    /* 	for (auto it = map_3D_2DW_set[i].first.begin(); it!= map_3D_2DW_set[i].first.end(); it++){ */
    /* 	  double charge = std::get<0>(map_2D_wt_charge[*it]); */
    /* 	  double charge_err = std::get<1>(map_2D_wt_charge[*it]); */
    /* 	  int n_divide =  map_2DW_3D_set[*it].size(); */
    /* 	  double scaling = charge/charge_err/n_divide; */
	  
    /* 	  ave_wt += it->second * scaling; */
    /* 	  ave_wc += it->first * scaling; */
    /* 	  ave_ws += scaling; */
    /* 	  //	  std::cout << i << " " << charge << " " << charge_err << " " << n_divide << std::endl; */
    /* 	} */

    /* 	double sum = 0; */
    /* 	if (ave_us > 0 && map_3D_2DU_set[i].second == 0) */
    /* 	  sum += pow(cur_t-ave_ut/ave_us,2)+pow(cur_u-ave_uc/ave_us,2); */
    /* 	if (ave_vs > 0 && map_3D_2DV_set[i].second == 0) */
    /* 	  sum += pow(cur_t-ave_vt/ave_vs,2)+pow(cur_v-ave_vc/ave_vs,2); */
    /* 	if (ave_ws > 0 && map_3D_2DW_set[i].second == 0) */
    /* 	  sum += pow(cur_t-ave_wt/ave_ws,2)+pow(cur_w-ave_wc/ave_ws,2); */
    /* 	//	 */
    /* 	if (sqrt(sum)>1.2) */
    /* 	  std::cout << cluster_id << " " << i << " " << sqrt(pow(cur_t-ave_ut/ave_us,2)+pow(cur_u-ave_uc/ave_us,2)) << " " << sqrt(pow(cur_t-ave_vt/ave_vs,2)+pow(cur_v-ave_vc/ave_vs,2)) << " " << sqrt(pow(cur_t-ave_wt/ave_ws,2)+pow(cur_w-ave_wc/ave_ws,2)) <<  " " << map_3D_2DU_set[i].second << " " << map_3D_2DV_set[i].second << " " << map_3D_2DW_set[i].second << std::endl; */
    /* 	//	  std::cout << cluster_id << " " << i << " " << sqrt(sum) << std::endl; */
	
    /*   } */
    /*   //   } */
    
    
    // std::cout << p.x << " " << p.y << " " << p.z << " " << path_wcps_vec.at(i).x << " " << path_wcps_vec.at(i).y << " " << path_wcps_vec.at(i).z << std::endl;
    /* if (cluster_id==18){ */
    /*   std::cout << slope_x * fine_tracking_path.at(6).x + offset_t << " "  */
    /* 		<< slope_yu *fine_tracking_path.at(6).y + slope_zu * fine_tracking_path.at(6).z + offset_u << " " */
    /* 		<< slope_yv *fine_tracking_path.at(6).y + slope_zv * fine_tracking_path.at(6).z + offset_v + 2400 << " " */
    /* 		<< slope_yw *fine_tracking_path.at(6).y + slope_zw * fine_tracking_path.at(6).z + offset_w + 4800<< std::endl; */
    /* } */
  }
  
}



void PR3DCluster::organize_ps_path(PointVector& ps_vec,  double low_dis_limit, std::vector<int>& record_vec){

  ps_vec.clear();
  // deal with the beginning ... 
  double dis1 = sqrt(pow(fine_tracking_path.at(1).x-fine_tracking_path.at(0).x,2)+pow(fine_tracking_path.at(1).y-fine_tracking_path.at(0).y,2)+pow(fine_tracking_path.at(1).z-fine_tracking_path.at(0).z,2));
  Point p1 = fine_tracking_path.at(0);
  if (dis1 > 0.3*units::cm){
    p1.x += (fine_tracking_path.at(0).x-fine_tracking_path.at(1).x)/dis1*(dis1-0.3*units::cm);
    p1.y += (fine_tracking_path.at(0).y-fine_tracking_path.at(1).y)/dis1*(dis1-0.3*units::cm);
    p1.z += (fine_tracking_path.at(0).z-fine_tracking_path.at(1).z)/dis1*(dis1-0.3*units::cm);
  }
  ps_vec.push_back(p1);
  record_vec.push_back(1);
  
  
  for (int i=0;i+1!=fine_tracking_path.size();i++){
    double dis = sqrt(pow(fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x,2)+pow(fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y,2)+pow(fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z,2));
    
    record_vec.push_back(0);
    if (dis >  low_dis_limit * 3.2){
      int npoints = std::round(dis/low_dis_limit/2.)-1;
      // double length = sqrt(pow(fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x,2)+pow(fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y,2)+pow(fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z,2));
      
      for (int j=0;j!=npoints;j++){
	Point p(fine_tracking_path.at(i).x + (fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x)/(npoints+1.)*(j+1.),
		fine_tracking_path.at(i).y + (fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y)/(npoints+1.)*(j+1.),
		fine_tracking_path.at(i).z + (fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z)/(npoints+1.)*(j+1.));
	ps_vec.push_back(p);
	record_vec.push_back(1);
      }
    }
  }
  record_vec.push_back(0);
  
  // deal with the end points ...

  double dis2 = sqrt(pow(fine_tracking_path.at(fine_tracking_path.size()-1).x-fine_tracking_path.at(fine_tracking_path.size()-2).x,2)+pow(fine_tracking_path.at(fine_tracking_path.size()-1).y-fine_tracking_path.at(fine_tracking_path.size()-2).y,2)+pow(fine_tracking_path.at(fine_tracking_path.size()-1).z-fine_tracking_path.at(fine_tracking_path.size()-2).z,2));
  Point p2 = fine_tracking_path.at(fine_tracking_path.size()-1);
  if (dis2 > 0.3*units::cm){
    p2.x += (fine_tracking_path.at(fine_tracking_path.size()-1).x-fine_tracking_path.at(fine_tracking_path.size()-2).x)/dis2*(dis2-0.3*units::cm);
    p2.y += (fine_tracking_path.at(fine_tracking_path.size()-1).y-fine_tracking_path.at(fine_tracking_path.size()-2).y)/dis2*(dis2-0.3*units::cm);
    p2.z += (fine_tracking_path.at(fine_tracking_path.size()-1).z-fine_tracking_path.at(fine_tracking_path.size()-2).z)/dis2*(dis2-0.3*units::cm);
  }
  ps_vec.push_back(p2);
  record_vec.push_back(1);
  
 
  
  fine_tracking_path.clear();
}

void PR3DCluster::merge_path(PointVector& fine_tracking_path_1st, PointVector& fine_tracking_path_2nd, std::vector<int>& record_vec, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DU_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DV_set, std::map<int,std::pair<std::set<std::pair<int,int>>, int> >& map_3D_2DW_set){

  fine_tracking_path.clear();
  int counter_1st = 0;
  int counter_2nd = 0;
  for (auto it = record_vec.begin(); it!=record_vec.end(); it++){
    if (*it==0){
      fine_tracking_path.push_back(fine_tracking_path_1st.at(counter_1st));
      counter_1st ++;
    }else if (*it==1){
      int sum = 0;
      
      if (map_3D_2DU_set[counter_2nd].first.size()!=0 && map_3D_2DU_set[counter_2nd].second == 0 ) sum ++;
      if (map_3D_2DV_set[counter_2nd].first.size()!=0 && map_3D_2DV_set[counter_2nd].second == 0 ) sum ++;
      if (map_3D_2DW_set[counter_2nd].first.size()!=0 && map_3D_2DW_set[counter_2nd].second == 0 ) sum ++;
      if ((counter_2nd == 0 || counter_2nd == fine_tracking_path_2nd.size()-1)){
	if (sum>1)
	  fine_tracking_path.push_back(fine_tracking_path_2nd.at(counter_2nd));
      }else{
	if (sum>0){
	  fine_tracking_path.push_back(fine_tracking_path_2nd.at(counter_2nd));	  
	}
      }
      counter_2nd ++;
    }

  }
}



void PR3DCluster::fine_tracking(std::map<int,std::map<const GeomWire*, SMGCSelection > >& global_wc_map, int num_pts_cut){

  // first round tracking with graph-based solution
  // if (path_wcps.size() < num_pts_cut) return;

  // organize trajectory point based on graph-based solution 
  double low_dis_limit = 0.6*units::cm;
  std::vector<WCPointCloud<double>::WCPoint> path_wcps_vec;
  organize_wcps_path(path_wcps_vec,low_dis_limit);

  //  std::cout << "temp " << std::endl;

  /* for (size_t i=0;i!=path_wcps_vec.size();i++){ */
  /*   std::cout << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << std::endl; */
  /* } */
  
  if (path_wcps_vec.size()>=num_pts_cut){
    // fill initial charge ...
    // ignore the dead channels here ... 
    std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_ut_charge;
    std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_vt_charge;
    std::map<std::pair<int,int>,std::tuple<double,double, int> > map_2D_wt_charge;

    //fill_2d_charge(global_wc_map, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);
    {
      std::vector<int> proj_channel;
      std::vector<int> proj_timeslice;
      std::vector<int> proj_charge;
      std::vector<int> proj_charge_err;
      std::vector<int> proj_flag;
      get_projection(proj_channel,proj_timeslice,proj_charge, proj_charge_err, proj_flag, global_wc_map);

      // std::cout << proj_charge.size() << " " << proj_flag.size() << std::endl;
      
      for (size_t i=0;i!=proj_channel.size();i++){
    	if (proj_channel.at(i)<2400){
    	  map_2D_ut_charge[std::make_pair(proj_channel.at(i),proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
    	}else if (proj_channel.at(i)<4800){
    	  map_2D_vt_charge[std::make_pair(proj_channel.at(i)-2400,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
    	}else{
    	  map_2D_wt_charge[std::make_pair(proj_channel.at(i)-4800,proj_timeslice.at(i))] = std::make_tuple(proj_charge.at(i),proj_charge_err.at(i), proj_flag.at(i));
    	}
      }
    }
    
    //    std::cout << "temp1 " << std::endl;
    
    // map index ... 
    // map 3D index to set of 2D points
    std::map<int,std::pair<std::set<std::pair<int,int>>, int> > map_3D_2DU_set;
    std::map<int,std::pair<std::set<std::pair<int,int>>, int> > map_3D_2DV_set;
    std::map<int,std::pair<std::set<std::pair<int,int>>, int> > map_3D_2DW_set;
    // map 2D points to 3D indices
    std::map<std::pair<int,int>,std::set<int>> map_2DU_3D_set;
    std::map<std::pair<int,int>,std::set<int>> map_2DV_3D_set;
    std::map<std::pair<int,int>,std::set<int>> map_2DW_3D_set;
    
    form_map_graph_based(path_wcps_vec, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
			 map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
			 map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set);


    
    /* for (int i=0;i!=path_wcps_vec.size();i++){ */
    /*   Point p(path_wcps_vec.at(i).x, path_wcps_vec.at(i).y, path_wcps_vec.at(i).z); */
    /*   std::cout << i <<  " " << map_3D_2DU_set[i].first.size() << " " << map_3D_2DV_set[i].first.size() << " " << map_3D_2DW_set[i].first.size() << " " << p << " " << path_wcps_vec.at(i).mcell << std::endl; */
    /* } */
    
    //    std::cout << "temp2 " << " " << path_wcps_vec.size() << std::endl;
    
    PointVector ps_vec;
    for (size_t i=0;i!=path_wcps_vec.size();i++){
      Point p(path_wcps_vec.at(i).x, path_wcps_vec.at(i).y, path_wcps_vec.at(i).z);
      ps_vec.push_back(p);
    }
    //  TPCParams& mp = Singleton<TPCParams>::Instance();
    //  double time_slice_width = mp.get_ts_width();
    //  double first_t_dis = path_wcps_vec.at(0).mcell->GetTimeSlice()*time_slice_width - path_wcps_vec.at(0).x;
    
    // trajectory fitting
    trajectory_fit(ps_vec,
    		   map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
    		   map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
    		   map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

    /* for (size_t i=0;i!=ps_vec.size();i++){ */
    /*   std::cout << i << " " << ps_vec.at(i).x/units::cm << " " << ps_vec.at(i).y/units::cm */
    /* 		<< " " << ps_vec.at(i).z/units::cm << std::endl; */
    /* } */

    //return;
    //    std::cout << "temp3 " << std::endl;
    
    /* for (int i=0;i+1!=path_wcps_vec.size();i++){ */
    /*   std::cout << sqrt(pow(path_wcps_vec.at(i+1).x-path_wcps_vec.at(i).x,2)+pow(path_wcps_vec.at(i+1).y-path_wcps_vec.at(i).y,2)+pow(path_wcps_vec.at(i+1).z-path_wcps_vec.at(i).z,2))/units::cm << " " << sqrt(pow(fine_tracking_path.at(i+1).x-fine_tracking_path.at(i).x,2)+pow(fine_tracking_path.at(i+1).y-fine_tracking_path.at(i).y,2)+pow(fine_tracking_path.at(i+1).z-fine_tracking_path.at(i).z,2))/units::cm << std::endl; */
    /* } */

    // second round fit ...
    path_wcps_vec.clear();
    for (size_t i=0;i!=fine_tracking_path.size();i++){
      /* std::cout << i << " " << fine_tracking_path.at(i).x/units::cm << */
      /* 	" " << fine_tracking_path.at(i).y/units::cm << " " */
      /* 		<< fine_tracking_path.at(i).z/units::cm << std::endl; */
      WCP::WCPointCloud<double>::WCPoint& wcp = point_cloud->get_closest_wcpoint(fine_tracking_path.at(i));
      if (wcp.index!= path_wcps_vec.back().index)
    	path_wcps_vec.push_back(wcp);
    }
    organize_wcps_vec_path(path_wcps_vec,low_dis_limit);

    //    std::cout << "temp4 " << " " << path_wcps_vec.size() << std::endl;

    /* for (size_t i=0;i!=path_wcps_vec.size();i++){ */
    /*   std::cout << path_wcps_vec.at(i).x/units::cm << " " */
    /* 		<< path_wcps_vec.at(i).y/units::cm << " " */
    /* 		<< path_wcps_vec.at(i).z/units::cm << std::endl; */
    /* } */
    
    if (path_wcps_vec.size()>1){
      form_map_graph_based(path_wcps_vec, map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge,
    			   map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
    			   map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set);


      
      ps_vec.clear();
      for (size_t i=0;i!=path_wcps_vec.size();i++){
    	Point p(path_wcps_vec.at(i).x, path_wcps_vec.at(i).y, path_wcps_vec.at(i).z);
    	ps_vec.push_back(p);
      }
      //      std::cout << "temp4_0" << " " << ps_vec.size() << std::endl;
      trajectory_fit(ps_vec,
    		     map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
    		     map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
    		     map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

      //      std::cout << "temp4_1" << " " << ps_vec.size() << std::endl;

      //      for (size_t i=0;i!=ps_vec.size();i++){
      //	std::cout << i << " " << ps_vec.at(i).x/units::cm << " " <<
      //	  ps_vec.at(i).y/units::cm << " " << ps_vec.at(i).z/units::cm << std::endl;
      //}
      
      PointVector fine_tracking_path_1st = fine_tracking_path;
      std::vector<int> record_vec;
      organize_ps_path(ps_vec,  low_dis_limit, record_vec);
      //std::cout << "temp4_1.3" << " " << ps_vec.size() << std::endl;
      //form association ...
      form_map_projection_based(ps_vec, map_3D_2DU_set, map_3D_2DV_set,  map_3D_2DW_set,
      				map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
      				map_2D_ut_charge,  map_2D_vt_charge, map_2D_wt_charge,
      				0.25, 0.9, 3, 4);
      
      //      std::cout << "temp4_1.5" << " " << ps_vec.size() << std::endl;
      // fit again ...
      trajectory_fit(ps_vec,
      		     map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set,
      		     map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set,
      		     map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge);

      //std::cout << "temp4_2" << std::endl;
      
      PointVector fine_tracking_path_2nd = fine_tracking_path;

      //      std::cout << fine_tracking_path_1st.size() << " " << fine_tracking_path_2nd.size() << " " << record_vec.size() << std::endl;
      
      merge_path(fine_tracking_path_1st, fine_tracking_path_2nd, record_vec,
    		 map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set);

      //      std::cout << "temp4_3" << std::endl;
    }
    
    
    /* // do a 2nd round iteration on the fit ... */
    /* ps_vec = fine_tracking_path; */
    /* // form association ... */
    /* form_map_projection_based(ps_vec, map_3D_2DU_set, map_3D_2DV_set,  map_3D_2DW_set, */
    /* 			      map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set, */
    /* 			      map_2D_ut_charge,  map_2D_vt_charge, map_2D_wt_charge, */
    /* 			      0.9, 0.9, 4, 4); */
    /* // fit again ... */
    /* trajectory_fit(ps_vec, */
    /* 		   map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set, */
    /* 		   map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set, */
    /* 		   map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge); */
    /* //Add a new round of fit to fill in the gap ... */
    /* PointVector fine_tracking_path_1st = fine_tracking_path; */
    /* std::vector<int> record_vec; */
    /* organize_ps_path(ps_vec,  low_dis_limit, record_vec); */
    /* //form association ... */
    /* form_map_projection_based(ps_vec, map_3D_2DU_set, map_3D_2DV_set,  map_3D_2DW_set, */
    /* 			      map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set, */
    /* 			      map_2D_ut_charge,  map_2D_vt_charge, map_2D_wt_charge, */
    /* 			      0.25, 0.9, 4, 4); */
    /* // fit again ... */
    /* trajectory_fit(ps_vec, */
    /* 		   map_3D_2DU_set, map_3D_2DV_set, map_3D_2DW_set, */
    /* 		   map_2DU_3D_set, map_2DV_3D_set, map_2DW_3D_set, */
    /* 		   map_2D_ut_charge, map_2D_vt_charge, map_2D_wt_charge); */
    

    /* PointVector fine_tracking_path_2nd = fine_tracking_path; */
    /* merge_path(fine_tracking_path_1st, fine_tracking_path_2nd, record_vec); */
    
    //   std::cout << fine_tracking_path_1st.size() << " " << fine_tracking_path_2nd.size() << " " << fine_tracking_path.size() << " " << record_vec.size() << std::endl;
    
  }else if (path_wcps_vec.size()>0){
    // for now, very short track, just copy .... 
    fine_tracking_path.clear();
    for (size_t i=0;i!=path_wcps_vec.size();i++){
      Point p;
      p.x = path_wcps_vec.at(i).x;
      p.y = path_wcps_vec.at(i).y;
      p.z = path_wcps_vec.at(i).z;
      fine_tracking_path.push_back(p);
    }
  }

  std::cout << "temp5 " << std::endl;
  
  // examine ...
  examine_path(low_dis_limit);
}


void PR3DCluster::examine_path(double low_dis_limit){
  PointVector temp_fine_tracking_path = fine_tracking_path;
  fine_tracking_path.clear();
  for (size_t i=0;i!=temp_fine_tracking_path.size();i++){
    if (i==0){
      fine_tracking_path.push_back(temp_fine_tracking_path.at(i));
    }else if (i+1==temp_fine_tracking_path.size()){
      double dis = sqrt(pow(temp_fine_tracking_path.at(i).x - fine_tracking_path.back().x,2)
  			+pow(temp_fine_tracking_path.at(i).y - fine_tracking_path.back().y,2)
  			+pow(temp_fine_tracking_path.at(i).z - fine_tracking_path.back().z,2));
      if (dis < low_dis_limit * 0.75){
	fine_tracking_path.pop_back();
	fine_tracking_path.push_back(temp_fine_tracking_path.at(i));
      }else{
	fine_tracking_path.push_back(temp_fine_tracking_path.at(i));
      }
    }else {
      double dis = sqrt(pow(temp_fine_tracking_path.at(i).x - fine_tracking_path.back().x,2)
  			+pow(temp_fine_tracking_path.at(i).y - fine_tracking_path.back().y,2)
  			+pow(temp_fine_tracking_path.at(i).z - fine_tracking_path.back().z,2));
      if (dis > low_dis_limit * 0.75){
	TVector3 v1(temp_fine_tracking_path.at(i).x-fine_tracking_path.back().x,
		    temp_fine_tracking_path.at(i).y-fine_tracking_path.back().y,
		    temp_fine_tracking_path.at(i).z-fine_tracking_path.back().z);
	TVector3 v2;
	double dis1 =  sqrt(pow(temp_fine_tracking_path.at(i).x - temp_fine_tracking_path.at(i+1).x,2)
			    +pow(temp_fine_tracking_path.at(i).y - temp_fine_tracking_path.at(i+1).y,2)
			    +pow(temp_fine_tracking_path.at(i).z - temp_fine_tracking_path.at(i+1).z,2));
	if (dis1 > low_dis_limit*0.75 || i+2==temp_fine_tracking_path.size()){
	  v2.SetXYZ(temp_fine_tracking_path.at(i).x-temp_fine_tracking_path.at(i+1).x,
		    temp_fine_tracking_path.at(i).y-temp_fine_tracking_path.at(i+1).y,
		    temp_fine_tracking_path.at(i).z-temp_fine_tracking_path.at(i+1).z);
	}else{
	  v2.SetXYZ(temp_fine_tracking_path.at(i).x-temp_fine_tracking_path.at(i+2).x,
		    temp_fine_tracking_path.at(i).y-temp_fine_tracking_path.at(i+2).y,
		    temp_fine_tracking_path.at(i).z-temp_fine_tracking_path.at(i+2).z);
	}
	fine_tracking_path.push_back(temp_fine_tracking_path.at(i));
      }
    }
  }

  temp_fine_tracking_path = fine_tracking_path;
  std::vector<double> temp_angle_vec;
  for (size_t i=0;i!=temp_fine_tracking_path.size();i++){
    if (i==0 || i+1 == temp_fine_tracking_path.size()){
      temp_angle_vec.push_back(3.1415926);
    }else{
      TVector3 v1(temp_fine_tracking_path.at(i).x - temp_fine_tracking_path.at(i-1).x, temp_fine_tracking_path.at(i).y - temp_fine_tracking_path.at(i-1).y, temp_fine_tracking_path.at(i).z - temp_fine_tracking_path.at(i-1).z);
      TVector3 v2(temp_fine_tracking_path.at(i).x - temp_fine_tracking_path.at(i+1).x, temp_fine_tracking_path.at(i).y - temp_fine_tracking_path.at(i+1).y, temp_fine_tracking_path.at(i).z - temp_fine_tracking_path.at(i+1).z);
      temp_angle_vec.push_back(v1.Angle(v2));
    }
  }

  fine_tracking_path.clear();
  for (size_t i=0;i!=temp_fine_tracking_path.size();i++){
    // if (temp_angle_vec.at(i)/3.1415926*180.>= 60)
      fine_tracking_path.push_back(temp_fine_tracking_path.at(i));
  }
  
  // std::cout << path_wcps_vec.size() << " " << fine_tracking_path.size() << std::endl;
}
