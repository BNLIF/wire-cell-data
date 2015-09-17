#include "WireCellData/WCVertex.h"


#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
//#include "Minuit2/MnMinos.h"
//#include "Minuit2/MnContours.h"
//#include "Minuit2/MnPlot.h"
#include "TVector3.h"

using namespace WireCell;

WCVertex::WCVertex(MergeSpaceCell& msc)
  : msc(&msc)
{
  center = msc.Get_Center();
}

WCVertex::~WCVertex(){
}

void WCVertex::reset_center(){
  center = msc->Get_Center();
}

void WCVertex::set_ky(WCTrack *track, double ky){
  if (tracks_ky.size() == 0 ){
    tracks_ky.resize(tracks.size());
  }
  while(tracks_ky.size() < tracks.size()){
    tracks_ky.push_back(0);
  }
  auto it = find(tracks.begin(),tracks.end(),track);
  tracks_ky.at(it-tracks.begin()) = ky;
}

void WCVertex::set_kz(WCTrack *track, double kz){
  if (tracks_kz.size() == 0 ){
    tracks_kz.resize(tracks.size());
  }
  while(tracks_kz.size() < tracks.size()){
    tracks_kz.push_back(0);
  }
  auto it = find(tracks.begin(),tracks.end(),track);
  tracks_kz.at(it-tracks.begin()) = kz;
}

double WCVertex::get_ky(WCTrack *track){
  auto it = find(tracks.begin(),tracks.end(),track);
  return tracks_ky.at(it-tracks.begin());
}

double WCVertex::get_kz(WCTrack *track){
  auto it = find(tracks.begin(),tracks.end(),track);
  return tracks_kz.at(it-tracks.begin());
}


//   return 0;
// }

double MyFCN::get_chi2(const std::vector<double> & xx) const{
  double chi2 = 0;
  double charge = 0;

  int ntracks = vertex->get_tracks().size();
  double x,y,z;
  double ky[100],kz[100];
  
  if (ntracks==1){
    x = vertex->get_msc()->Get_Center().x;
    y = xx[0];
    z = xx[1];
    ky[0] = xx[2];
    kz[0] = xx[3];
  }else{
    x = xx[0];
    y = xx[1];
    z = xx[2];
    for (int i=0;i!=ntracks;i++){
      ky[i] = xx[3+2*i];
      kz[i] = xx[3+2*i+1];
    }
  }

  double x0 = vertex->get_msc()->Get_Center().x;

  for (int i=0;i!=ntracks;i++){
    WCTrack *track = vertex->get_tracks().at(i);
    MergeSpaceCellSelection cells = track->get_all_cells();

    double x01=0;
    int flag_check_range=1;
    //auto it = find(track->get_end_scells().begin(),track->get_end_scells().end(),mscell);
    if (vertex->get_msc() == track->get_end_scells().at(0)){
      x01 = track->get_end_scells().at(1)->Get_Center().x;
    }else if (vertex->get_msc() == track->get_end_scells().at(1)){
      x01 = track->get_end_scells().at(0)->Get_Center().x;
    }else{
      flag_check_range = 0;
    }

    //calculate medium
    std::vector<double> dy_all,dz_all;
    
    for (int j=0;j!=cells.size();j++){
      MergeSpaceCell *mscell = cells.at(j);
      MergeSpaceCell *prev_mscell = mscell;
      MergeSpaceCell *next_mscell = mscell;
      
      if (cells.size()>=3){
	if (j==0){
	  prev_mscell = cells.at(2);
	  next_mscell = cells.at(1);
	}else if (j==cells.size()-1){
	  prev_mscell = cells.at(cells.size()-2);
	  next_mscell = cells.at(cells.size()-3);
	}else{
	  prev_mscell = cells.at(j-1);
	  next_mscell = cells.at(j+1);
	}
      }


      double xc = mscell->Get_Center().x;

      int flag = 0;
      if (vertex->get_fit_type() == 1){
	if (fabs(xc-x0)/units::cm < 5)
	  flag = 1;
      }else if (vertex->get_fit_type()==2){
	if ((fabs(xc-x0)/units::cm < 5 && fabs(xc-x0)/units::cm > 0.5 && fit_flag == 0) || (fabs(xc-x0)/units::cm < 6 && fabs(xc-x0)/units::cm > 1.0 && fit_flag != 0))
	  flag = 1;
      }

      if (flag==1){
	double dy = sqrt((pow(mscell->get_dy(),2) + pow(prev_mscell->get_dy(),2) + pow(next_mscell->get_dy(),2))/3.);
	double dz = sqrt((pow(mscell->get_dz(),2) + pow(prev_mscell->get_dz(),2) + pow(next_mscell->get_dz(),2))/3.);
	
	// double dy = mscell->get_dy();
	// double dz = mscell->get_dz();
	
	
	if (dy == 0) dy = 0.3 * units::cm/2.;
	if (dz == 0) dz = 0.3 * units::cm/2.;
	
	dy_all.push_back(dy);
	dz_all.push_back(dz);
	
      }

    }
    
    double dy_ave, dz_ave;
    size_t n_ave = dy_all.size()/2.;
    //   std::cout << n_ave << " " << dy_all.size() << std::endl;
    if (dy_all.size() > 0 ){
      nth_element(dy_all.begin(),dy_all.begin()+n_ave,dy_all.end());
      dy_ave = dy_all[n_ave];
      nth_element(dz_all.begin(),dz_all.begin()+n_ave,dz_all.end());
      dz_ave = dz_all[n_ave];
    }else{
      dy_ave = 0.3 * units::cm/2.;
      dz_ave = 0.3 * units::cm/2.;

    }

    for (int j=0;j!=cells.size();j++){
      MergeSpaceCell *mscell = cells.at(j);
      MergeSpaceCell *prev_mscell = mscell;
      MergeSpaceCell *next_mscell = mscell;
      
      if (cells.size()>=3){
	if (j==0){
	  prev_mscell = cells.at(2);
	  next_mscell = cells.at(1);
	}else if (j==cells.size()-1){
	  prev_mscell = cells.at(cells.size()-2);
	  next_mscell = cells.at(cells.size()-3);
	}else{
	  prev_mscell = cells.at(j-1);
	  next_mscell = cells.at(j+1);
	}
      }


      double xc = mscell->Get_Center().x;

      


      
      // if (vertex->get_fit_type() == 1){
      // 	if (fabs(xc-x0)/units::cm < 5)
      // 	  flag = 1;
      // }else if (vertex->get_fit_type()==2){
      // 	if ((fabs(xc-x0)/units::cm < 5 && fabs(xc-x0)/units::cm > 0.5 && fit_flag == 0) || (fabs(xc-x0)/units::cm < 6 && fabs(xc-x0)/units::cm > 1.0 && fit_flag != 0))
      // 	  flag = 1;
      // }
      


      // if ( (fabs(xc-x0)/units::cm < 5 && ntracks < 2) 
      // 	   || (fabs(xc-x0)/units::cm < 5 && fabs(xc-x0)/units::cm > 0.5 && ntracks >=2)){
      // 	flag = 1;
      // }
      // if (ntracks==1){
      // 	if (fabs(xc-x0)/units::cm > 0.1 && )
      // 	  flag = 1;
      // }else if (ntracks>1){
      // 	if (fabs(xc-x0)/units::cm > 0.7 && fabs(xc-x0)/units::cm < 5)
      // 	  flag = 1;
      // }
      // if (it != track->get_end_scells().end()){
      // 	flag = 0;
      // }

      
      
      double x1 = mscell->Get_Center().x;
      double y1 = mscell->Get_Center().y;
      double z1 = mscell->Get_Center().z;
      //double q = mscell->Get_Charge();
      double q = 1;
      
      double dy = sqrt((pow(mscell->get_dy(),2) + pow(prev_mscell->get_dy(),2) + pow(next_mscell->get_dy(),2))/3.);
      double dz = sqrt((pow(mscell->get_dz(),2) + pow(prev_mscell->get_dz(),2) + pow(next_mscell->get_dz(),2))/3.);
      
      // double dy = mscell->get_dy();
      // double dz = mscell->get_dz();
      
      
      if (dy == 0) dy = 0.3 * units::cm/2.;
      if (dz == 0) dz = 0.3 * units::cm/2.;
      

      int flag = 0;
      if (vertex->get_fit_type() == 1){
      	if (fabs(xc-x0)/units::cm < 5)
      	  flag = 1;
      }else if (vertex->get_fit_type()==2){
      	if (((fabs(xc-x0)/units::cm < 5 && fit_flag == 0) || (fabs(xc-x0)/units::cm < 6 && fit_flag != 0))&&(dy<dy_ave*2&&dz<dz_ave*2))
      	  flag = 1;
      }


      double a,b,c,d;
      
      a = y - ky[i]*x;
      b = ky[i];
      c = z - kz[i]*x;
      d = kz[i];
      
      double x2 = (-a*b-c*d+b*y1+d*z1+x1)/(1+b*b+d*d);
      double y2 = a + b*x2;
      double z2 = c + d*x2;
      
      //	chi2 += (pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2))*q/0.15/0.15*3./units::cm/units::cm;
      if (flag ==1 ){	
	chi2 += (pow(x2-x1,2)/0.16/0.16*3./units::cm/units::cm +  
		 pow(y2-y1,2)/pow(dy,2)*3 +
		 pow(z2-z1,2)/pow(dz,2)*3
		 )*q;
	
	//  	std::cout << i << " " << j << " " << x << " " << y << " " << z << " " << ky[i] << " " << kz[i] << std::endl;
	
	charge += q;
	
	
	// std::cout << flag_check_range << " " << i << " " << x1 << " " << y1 << " " << z1 << " " 
	// 	  << x << " " << x0 << " " << x01 << std::endl;
	
	if (flag_check_range == 1){
	  if ( (x01-x0)*(x1-x)<0){
	    chi2 += (pow(x2-x1,2)/0.16/0.16*3./units::cm/units::cm +  
		     pow(y2-y1,2)/pow(dy,2)*3 +
		     pow(z2-z1,2)/pow(dz,2)*3
		     )*10;
	  }
	}
      }else{
	// if (fabs(x2-x1) <= 0.16*units::cm && 
	//     fabs(y2-y1) <= fabs(dy) &&
	//     fabs(z2-z1) <= fabs(dz) ){
	// }else{
	//   chi2 += q;
	// }
      }


      // if (flag == 1 ){
      // 	// 1.0-5 cm
	
      // 	//std::cout << xc/units::cm << " " << x0/units::cm << std::endl;

      // 	SpaceCellSelection scells = mscell->Get_all_spacecell();
      // 	for (int k=0;k!=scells.size();k++){
      // 	  double x1 = scells.at(k)->x();
      // 	  double y1 = scells.at(k)->y();
      // 	  double z1 = scells.at(k)->z();
      // 	  double q = scells.at(k)->q();

      // 	  double a,b,c,d;
	  
      // 	  a = y - ky[i]*x;
      // 	  b = ky[i];
      // 	  c = z - kz[i]*x;
      // 	  d = kz[i];

      // 	  double x2 = (-a*b-c*d+b*y1+d*z1+x1)/(1+b*b+d*d);
      // 	  double y2 = a + b*x2;
      // 	  double z2 = c + d*x2;

      // 	  // std::cout << x1/units::cm << " " << x2/units::cm << " " << 
      // 	  //   y1/units::cm << " " << y2/units::cm << " " << 
      // 	  //   z1/units::cm << " " << z2/units::cm << " " <<  q << 
      // 	  //   std::endl;
      // 	  // TVector3 n(1,ky[i],kz[i]);
      // 	  // TVector3 p(x-x1,y-x1,z-x1);
      // 	  // TVector3 v = p.Cross(n);
      // 	  // double dis = v.Mag()/n.Mag();
	  
      // 	  chi2 += (pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2))*q/0.15/0.15*3./units::cm/units::cm;
      // 	  charge += q;
      // 	}
      //}
    }
  }

  //std::cout << chi2 << " " << charge << std::endl;
  
  return chi2/charge;
}

double MyFCN::operator() (const std::vector<double> & xx) const{
  return get_chi2(xx);
  
} 


bool WCVertex::FindVertex(int flag){

  // hack ...
  // if (tracks.size()>3){
  //   //    tracks.erase(tracks.begin()+4);
  //   tracks.erase(tracks.begin()+3);
  //   tracks.erase(tracks.begin()+2);
  //   tracks.erase(tracks.begin()+1);
  //   tracks.erase(tracks.begin()+0);    
  // }

  WCTrackSelection bad_tracks;
  for (int i=0;i!=tracks.size();i++){
    if (tracks.at(i)->get_all_cells().size()==1)
      bad_tracks.push_back(tracks.at(i));
  }
  for (int i=0;i!=bad_tracks.size();i++){
    auto it = find(tracks.begin(),tracks.end(),bad_tracks.at(i));
    tracks.erase(it);
  }




  MyFCN fcn(this,flag);
  int ntracks = tracks.size();
  
  std::vector<double> variable;
  int npar;


  


  
  if (ntracks == 1){
    // 4 parameters, (2 for vertex y,z) and two directions
    npar = 4;
    variable.push_back(center.y);
    variable.push_back(center.z);
    
    WCTrack* track = tracks.at(0);
    MergeSpaceCell *a1 = track->get_end_scells().at(0);
    MergeSpaceCell *a2 = track->get_end_scells().at(1);

    // float dis1 = fabs(a1->Get_Center().x-center.x);
    // float dis2 = fabs(a2->Get_Center().x-center.x);
    float ky,kz;
    
    ky = (a2->Get_Center().y-a1->Get_Center().y)/(a2->Get_Center().x-a1->Get_Center().x);
    kz = (a2->Get_Center().z-a1->Get_Center().z)/(a2->Get_Center().x-a1->Get_Center().x);
    
    //std::cout << "abc: " << center.y/units::cm << " " << center.z/units::cm << " " << ky << " " << kz << std::endl; 
    
    variable.push_back(ky);
    variable.push_back(kz);
  }else{
    // 3 parameters for vertex and 2*n for track direction, each for a track
    npar = 3 + 2*ntracks;

    variable.push_back(center.x);
    variable.push_back(center.y);
    variable.push_back(center.z);
    
    for (int i=0;i!=ntracks;i++){
      WCTrack* track = tracks.at(i);
      MergeSpaceCell *a1 = track->get_end_scells().at(0);
      MergeSpaceCell *a2 = track->get_end_scells().at(1);
      
      
      float ky,kz;
      
      ky = (a2->Get_Center().y-a1->Get_Center().y)/(a2->Get_Center().x-a1->Get_Center().x);
      kz = (a2->Get_Center().z-a1->Get_Center().z)/(a2->Get_Center().x-a1->Get_Center().x);
     
      //std::cout << "abc: " << ky << " " << kz << std::endl; 
      
      variable.push_back(ky);
      variable.push_back(kz);
    }
  }
  
  //figure out fit_type;
  // fit 1, fit everything
  // fit 2 exclude vertex and fit
  
  fit_type = 2;
  int num_tracks = 0;
  for (int i=0;i!=tracks.size();i++){
    int num_cells = 0;
    for (int j=0;j!=tracks.at(i)->get_all_cells().size();j++){
      if (flag == 0){
      if (fabs(msc->Get_Center().x-tracks.at(i)->get_all_cells().at(j)->Get_Center().x)/units::cm > 0.5)
	num_cells ++;
      }else{
	if (fabs(msc->Get_Center().x-tracks.at(i)->get_all_cells().at(j)->Get_Center().x)/units::cm > 1.0)
	  num_cells ++;
      }
    }
    if (num_cells >= 4) num_tracks ++;
  }
  if (num_tracks<2) fit_type = 1;

  //std::cout << fit_type << " " << num_tracks << " " << ntracks << std::endl;
  

  //fit_type = 1;


  ROOT::Minuit2::MnUserParameters upar;
  for (int i=0;i!=npar;i++){
    upar.Add(Form("x_%d",i),variable.at(i),0.001);
  }
  // create MIGRAD minimizer
  ROOT::Minuit2::MnMigrad migrad(fcn, upar);

  // minimize
  ROOT::Minuit2::FunctionMinimum min = migrad();
  const ROOT::Minuit2::MnUserParameterState& state = min.UserState();
  
  tracks_ky.clear();
  tracks_kz.clear();
  if (ntracks == 1){
    center.y = state.Value(0);
    center.z = state.Value(1);
    tracks_ky.push_back(state.Value(2));
    tracks_kz.push_back(state.Value(3));
  }else{
    center.x = state.Value(0);
    center.y = state.Value(1);
    center.z = state.Value(2);

    

    for (int i=0;i!=ntracks;i++){
      tracks_ky.push_back(state.Value(3+2*i));
      tracks_kz.push_back(state.Value(3+2*i+1));

      //std::cout << state.Value(3+2*i) << " " << state.Value(3+2*i+1) << std::endl;
    }
  }
  
  //std::cout << min.IsValid() << std::endl;
  

  // output
  //  std::cout << fcn.get_chi2(variable) << " Minimum: " << min.Fval() << " " << state.Value(0)/units::cm << " " <<state.Value(1)/units::cm << " " << state.Value(2)/units::cm << std::endl;
  // std::cout << "abc1: " << min.Fval() << " " << fcn.get_chi2(variable) << " " << fcn.get_chi2(variable1) << " " << state.Value(0)/units::cm << " " <<state.Value(1)/units::cm  << " " << state.Value(2) << " " << state.Value(3) << std::endl;
  //std::cout << fcn.get_chi2(variable) << std::endl;


  fit_success = min.IsValid();
  fit_chi2 = min.Fval();
  

  //return 0;
  return min.IsValid();
}




Point WCVertex::Center(){
  
  return center;
}

void WCVertex::Add(WCTrack* track){
  tracks.push_back(track);
}

void WCVertex::OrganizeTracks(){
  
  if (tracks.size()>1){
    
    // find a track where the vertex is one of the end track
    WCTrack *end_track;
    WCTrackSelection end_tracks;
    for (int i = 0;i!=tracks.size();i++){
      auto it = find(tracks.at(i)->get_end_scells().begin(),
		     tracks.at(i)->get_end_scells().end(),
		     msc);
      if (it != tracks.at(i)->get_end_scells().end()){
	end_track = tracks.at(i);
	end_tracks.push_back(end_track);
      }
    }
    
    
    //    MergeClusterTrack& mct = end_track->get_mct();
    // find out where is the mct
    
    int flag;
    Point p0 = msc->Get_Center();
    Point p1 = end_track->get_all_cells().front()->Get_Center();
    Point p2 = end_track->get_all_cells().back()->Get_Center();
    float dis1 = sqrt(pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2));
    float dis2 = sqrt(pow(p0.x-p2.x,2)+pow(p0.y-p2.y,2)+pow(p0.z-p2.z,2));
    
    if (dis1 < dis2){
      flag = 1;
    }else{
      flag = -1;
    }
    
    int time_length = end_track->get_all_cells().size();
    MergeSpaceCell* nvertex = msc;
    
    if (flag ==1){
      int flag_common = 0;      
      for (int i=0;i!=time_length;i++){
	MergeSpaceCell *tcell = end_track->get_all_cells().at(i);
	if (this->CheckContain(tcell)){
	  flag_common = 1;
	  if (tcell->Get_Center().x != nvertex->Get_Center().x){
	    nvertex = tcell;
	  }else{
	    if (tcell->Get_all_spacecell().size() > nvertex->Get_all_spacecell().size()){
	      nvertex = tcell;
	    } 
	  }
	}
    	if (flag_common == 0)
    	  break;
      }
      
    }else{
      int flag_common = 0;      
      for (int i=0;i!=time_length;i++){
	MergeSpaceCell *tcell = end_track->get_all_cells().at(time_length-1-i);
	if (this->CheckContain(tcell)){
	  flag_common = 1;
	  if (tcell->Get_Center().x != nvertex->Get_Center().x){
	    nvertex = tcell;
	  }else{
	    if (tcell->Get_all_spacecell().size() > nvertex->Get_all_spacecell().size()){
	      nvertex = tcell;
	    } 
	  }
	}
    	if (flag_common == 0)
    	  break;
      }
    }

    
    for (int i=0;i!=end_tracks.size();i++){
      end_tracks.at(i)->ReplaceEndCell(msc,nvertex);
      end_tracks.at(i)->ModifyCells();
    }


    msc = nvertex;
    center = msc->Get_Center();
   
    // for (int i=0;i!=tracks.size();i++){
    //   tracks.at(i)->ModifyCells();
    // }
    
  }
  
}





WCTrackSelection WCVertex::BreakTracksAngle(WCTrackSelection& finished_tracks){
  WCTrackSelection result_tracks;

  for (int i=0;i!=tracks.size();i++){
    WCTrack *track = tracks.at(i);

    auto it = find(finished_tracks.begin(),finished_tracks.end(),track);
    if (it == finished_tracks.end()){
      finished_tracks.push_back(track);

      MergeSpaceCellSelection all_cells = track->get_all_cells();
      float dis1 = fabs(msc->Get_Center().x - all_cells.at(0)->Get_Center().x);
      float dis2 = fabs(msc->Get_Center().x - all_cells.at(all_cells.size()-1)->Get_Center().x);
      

      // try to find all the break points
      MergeSpaceCellSelection break_cells;
      //int prev_break_point = 0;

      MergeSpaceCell *start_cell = all_cells.front();
      MergeSpaceCell *end_cell = all_cells.back();
      MergeSpaceCell *prev_break_cell = all_cells.front();

      for (int j=0;j!=all_cells.size();j++){
	
	MergeSpaceCell* curr_cell;
	MergeSpaceCell* next_cell;
	MergeSpaceCell* prev_cell;

	MergeSpaceCellSelection curr_cells;
	MergeSpaceCellSelection prev_cells;
	MergeSpaceCellSelection next_cells;
	
	if (dis1 < dis2){
	  curr_cell = all_cells.at(j);
	}else{
	  curr_cell = all_cells.at(all_cells.size()-1-j);
	}
	
	//	curr_cells.push_back(curr_cell);
	for (int k=-5;k!=5;k++){
	  int num = j+k;
	  if (num >=0&&num <= all_cells.size()-1){
	    if (fabs(all_cells.at(num)->Get_Center().x - curr_cell->Get_Center().x)<0.1*units::cm){
	      curr_cells.push_back(all_cells.at(num));
	    }
	  }
	}

	
	// if (fabs(curr_cell->Get_Center().x/units::cm - 54.08)<2.0 )
	//   std::cout << j << " " << curr_cell->Get_Center().x/units::cm << std::endl;


	if (fabs(curr_cell->Get_Center().x - start_cell->Get_Center().x) >=5 * 0.31*units::cm 
	    && fabs(curr_cell->Get_Center().x - end_cell->Get_Center().x) >=5 * 0.31*units::cm 
	    && fabs(curr_cell->Get_Center().x - prev_break_cell->Get_Center().x) >=5 * 0.31*units::cm ){

	  // MergeSpaceCellSelection temp_msc;

	  if (dis1 < dis2){
	    // temp_msc.push_back(all_cells.at(j-1));
	    // temp_msc.push_back(all_cells.at(j-2));
	    // temp_msc.push_back(all_cells.at(j-3));
	    prev_cell = all_cells.at(j-4);
	    int k = 0;
	    while(fabs(prev_cell->Get_Center().x - curr_cell->Get_Center().x) < 4*0.31*units::cm){
	      k++;
	      if (j-4-k>=0&&j-4-k<all_cells.size()){
		prev_cell = all_cells.at(j-4-k);
	      }else{
		break;
	      }
	      // temp_msc.push_back(prev_cell);
	    }

	    for( int kk = -5;kk!=5;kk++){
	      int num = j-4-k+kk;
	      if (num >=0&&num <= all_cells.size()-1){
		if (fabs(all_cells.at(num)->Get_Center().x- prev_cell->Get_Center().x)<0.1*units::cm){
		  prev_cells.push_back(all_cells.at(num));
		}
	      }
	    }


	    // temp_msc.push_back(all_cells.at(j+1));
	    // temp_msc.push_back(all_cells.at(j+2));
	    // temp_msc.push_back(all_cells.at(j+3));
	    next_cell = all_cells.at(j+4);
	    
	    k = 0;
	    while(fabs(next_cell->Get_Center().x - curr_cell->Get_Center().x) < 4*0.31*units::cm){
	      k++;
	      if (j+4+k>=0 && j+4+k < all_cells.size()){
		next_cell = all_cells.at(j+4+k);
	      }else{
		break;
	      }
	      //temp_msc.push_back(next_cell);
	    }
	    
	    for( int kk = -5;kk!=5;kk++){
	      int num = j+4+k+kk;
	      if (num >=0&&num <= all_cells.size()-1){
		if (fabs(all_cells.at(num)->Get_Center().x - next_cell->Get_Center().x)<0.1*units::cm){
		  next_cells.push_back(all_cells.at(num));
		}
	      }
	    }


	  }else{
	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j+1));
	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j+2));
	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j+3));
	    prev_cell = all_cells.at(all_cells.size()-1-j+4);
	    int k=0;
	    while(fabs(prev_cell->Get_Center().x - curr_cell->Get_Center().x) < 4*0.31*units::cm){
	      k++;
	      if (all_cells.size()-1-j+4+k>=0&&all_cells.size()-1-j+4+k < all_cells.size()){
		prev_cell = all_cells.at(all_cells.size()-1-j+4+k);
	      }else{
		break;
	      }
	      //temp_msc.push_back(prev_cell);
	      
	    }
	    
	    for( int kk = -5;kk!=5;kk++){
	      int num = all_cells.size()-1-j+4+k+kk;
	      if (num >=0&&num <= all_cells.size()-1){
		if (fabs(all_cells.at(num)->Get_Center().x - prev_cell->Get_Center().x)<0.1*units::cm){
		  prev_cells.push_back(all_cells.at(num));
		}
	      }
	    }


	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j-1));
	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j-2));
	    // temp_msc.push_back(all_cells.at(all_cells.size()-1-j-3));
	    next_cell = all_cells.at(all_cells.size()-1-j-4);
	    k = 0;
	     while(fabs(next_cell->Get_Center().x - curr_cell->Get_Center().x) < 4*0.31*units::cm){
	      k++;
	      if (all_cells.size()-1-j-4-k>=0&&all_cells.size()-1-j-4-k<all_cells.size()){
		prev_cell = all_cells.at(all_cells.size()-1-j-4-k);
	      }else{
		break;
	      }
	      //temp_msc.push_back(next_cell);
	      
	     }

	     for( int kk = -5;kk!=5;kk++){
	      int num = all_cells.size()-1-j-4-k+kk;
	      if (num >=0&&num <= all_cells.size()-1){
		if (fabs(all_cells.at(num)->Get_Center().x - next_cell->Get_Center().x)<0.1*units::cm){
		  next_cells.push_back(all_cells.at(num));
		}
	      }
	    }

	  }




	  Point p1;
	  p1.x = 0; p1.y = 0; p1.z = 0;
	  int np1=0;

	  for (int kk=0;kk!=prev_cells.size();kk++){
	    p1.x += prev_cells.at(kk)->Get_Center().x * prev_cells.at(kk)->Get_all_spacecell().size();
	    p1.y += prev_cells.at(kk)->Get_Center().y * prev_cells.at(kk)->Get_all_spacecell().size();
	    p1.z += prev_cells.at(kk)->Get_Center().z * prev_cells.at(kk)->Get_all_spacecell().size();
	    np1 += prev_cells.at(kk)->Get_all_spacecell().size();
	  }
	  p1.x/=np1;
	  p1.y/=np1;
	  p1.z/=np1;

	  
	  Point p2;
	  p2.x = 0; p2.y = 0; p2.z = 0;
	  int np2=0;

	  for (int kk=0;kk!=curr_cells.size();kk++){
	    p2.x += curr_cells.at(kk)->Get_Center().x * curr_cells.at(kk)->Get_all_spacecell().size();
	    p2.y += curr_cells.at(kk)->Get_Center().y * curr_cells.at(kk)->Get_all_spacecell().size();
	    p2.z += curr_cells.at(kk)->Get_Center().z * curr_cells.at(kk)->Get_all_spacecell().size();
	    np2 += curr_cells.at(kk)->Get_all_spacecell().size();
	  }
	  p2.x/=np2;
	  p2.y/=np2;
	  p2.z/=np2;


	  Point p3;
	  p3.x = 0; p3.y = 0; p3.z = 0;
	  int np3=0;

	  for (int kk=0;kk!=next_cells.size();kk++){
	    p3.x += next_cells.at(kk)->Get_Center().x * next_cells.at(kk)->Get_all_spacecell().size();
	    p3.y += next_cells.at(kk)->Get_Center().y * next_cells.at(kk)->Get_all_spacecell().size();
	    p3.z += next_cells.at(kk)->Get_Center().z * next_cells.at(kk)->Get_all_spacecell().size();
	    np3 += next_cells.at(kk)->Get_all_spacecell().size();
	  }
	  p3.x/=np3;
	  p3.y/=np3;
	  p3.z/=np3;

	  

	  // Point p2 = curr_cell->Get_Center();
	  // Point p3 = next_cell->Get_Center();
	  
	  double ky = (p3.y-p1.y)/(p3.x-p1.x);
	  double kz = (p3.z-p1.z)/(p3.x-p1.x);

	  double yp = ky * (p2.x-p1.x) + p1.y;
	  double zp = kz * (p2.x-p1.x) + p1.z;



	  double dy = sqrt((pow(curr_cell->get_dy(),2)+pow(prev_cell->get_dy(),2)+pow(next_cell->get_dy(),2))/3.);
	  double dz = sqrt((pow(curr_cell->get_dz(),2)+pow(prev_cell->get_dz(),2)+pow(next_cell->get_dz(),2))/3.);
	  
	  if (dy < curr_cell->get_dy()) dy = curr_cell->get_dy();
	  if (dz < curr_cell->get_dz()) dz = curr_cell->get_dz();
	  if (dy == 0) dy = 0.3 * units::cm;
	  if (dz == 0) dz = 0.3 * units::cm;
	  dy = sqrt(dy*dy + pow(0.3*units::cm,2));
	  dz = sqrt(dz*dz + pow(0.3*units::cm,2));

	  if (curr_cells.size()>1){
	    double ymax = -1e9;
	    double ymin = 1e9;
	    double zmax = -1e9;
	    double zmin = 1e9;
	    for (int qx = 0;qx!=curr_cells.size();qx++){
	      if (curr_cells.at(qx)->get_maxy() > ymax) ymax = curr_cells.at(qx)->get_maxy();
	      if (curr_cells.at(qx)->get_maxz() > zmax) zmax = curr_cells.at(qx)->get_maxz();
	      if (curr_cells.at(qx)->get_miny() < ymin) ymin = curr_cells.at(qx)->get_miny();
	      if (curr_cells.at(qx)->get_minz() < zmin) zmin = curr_cells.at(qx)->get_minz();
	    }
	    
	    if (dy < fabs(ymax - ymin)/2.) dy = fabs(ymax - ymin)/2.;
	    if (dz < fabs(zmax - zmin)/2.) dz = fabs(zmax - zmin)/2.;

	  }
	  
	  

	  
	  





	  double dis_sigma2 = pow(yp-p2.y,2)/pow(dy,2)+pow(zp-p2.z,2)/pow(dz,2);
	  
	  float theta1_old = atan((p2.y-p1.y)/(p2.x-p1.x))/3.1415926*180.;
	  float theta2_old = atan((p2.z-p1.z)/(p2.x-p1.x))/3.1415926*180.;

	  float theta1_new = atan((p3.y-p2.y)/(p3.x-p2.x))/3.1415926*180.;
	  float theta2_new = atan((p3.z-p2.z)/(p3.x-p2.x))/3.1415926*180.;

	  if (sqrt(pow(theta1_new-theta1_old,2)+pow(theta2_new-theta2_old,2))>30. && dis_sigma2 > 1.4
	      && 3*curr_cell->get_dy() > prev_cell->get_dy() 
	      && 3*curr_cell->get_dy() > next_cell->get_dy() 
	      && 3*curr_cell->get_dz() > prev_cell->get_dz() 
	      && 3*curr_cell->get_dz() > next_cell->get_dz() 
	      ){
	    // std::cout << curr_cell->get_dy() << " " << prev_cell->get_dy()  << " " << next_cell->get_dy()
	    // 	      << " " << curr_cell->get_dz() << " " << prev_cell->get_dz()  << " " << next_cell->get_dz()
	    // 	      << std::endl;


	    prev_break_cell = curr_cell;
	    break_cells.push_back(curr_cell);

	    std::cout << "Break Angle: " << sqrt(pow(theta1_new-theta1_old,2)+pow(theta2_new-theta2_old,2)) << " " << dis_sigma2 << " " << 
	      curr_cell->Get_Center().x/units::cm << " " << curr_cell->Get_Center().y/units::cm << 
	      " " << curr_cell->Get_Center().z/units::cm << " " <<  curr_cell->get_dy() << " " 
		      << curr_cell->get_dz() << " "  << 
	      prev_cells.size() << " " << curr_cells.size() << " " << next_cells.size() << " " << 
	      // all_cells.at(0)->Get_Center().x/units::cm << " " << all_cells.back()->Get_Center().x/units::cm << " " << 
	      // prev_cell->Get_Center().x/units::cm << " " << prev_cell->Get_Center().y/units::cm << 
	      // " " << prev_cell->Get_Center().z/units::cm << " " <<  prev_cell->get_dy() << " " 
	      // 	      << prev_cell->get_dz() << " "  << 
	      // next_cell->Get_Center().x/units::cm << " " << next_cell->Get_Center().y/units::cm << 
	      // " " << next_cell->Get_Center().z/units::cm << " " <<  next_cell->get_dy() << " " 
	      // 	      << next_cell->get_dz() << " "  << 
	      // (yp - p2.y)/units::cm << " " << (zp-p2.z)/units::cm << " " << dy/units::cm << " " << dz/units::cm << " " << pow(yp-p2.y,2)/pow(dy,2)+pow(zp-p2.z,2)/pow(dz,2) << " " << 
	      std::endl;
	  }


	  // for (int k=0;k!= temp_msc.size();k++){
	  //   Point p = temp_msc.at(k)->Get_Center();
	  //   double yp1 = ky * (p.x-p1.x) + p1.y;
	  //   double zp1 = kz * (p.x-p1.x) + p1.z;
	    
	  //   double_t dis_sigma1 = sqrt(pow(yp-p.y,2)/pow(dy,2)+pow(zp-p.z,2)/pow(dz,2));
	  // }
	  
	  

	}
      }

     
      if (break_cells.size() > 0){
	track->get_end_scells().clear();
	track->get_all_cells().clear();
	
	//      std::cout << all_cells.size() << std::endl;
	WCTrack *track1 = track;
	result_tracks.push_back(track1);
	
	// use break points to form tracks
	for (int j=0;j!=all_cells.size();j++){
	  MergeSpaceCell *curr_cell;
	  if (dis1 < dis2){
	    curr_cell = all_cells.at(j);
	  }else{
	    curr_cell = all_cells.at(all_cells.size()-1-j);
	  }
	  
	  // push any way
	  if (j==0 || j == all_cells.size()-1){
	    track1->get_end_scells().push_back(curr_cell);
	  }
	  
	  track1->get_all_cells().push_back(curr_cell);
	  auto it = find(break_cells.begin(),break_cells.end(),curr_cell);
	  if (it != break_cells.end() ){
	    track1->get_end_scells().push_back(curr_cell);
	    track1 = new WCTrack(track->get_mct());
	    track1->get_end_scells().clear();
	    track1->get_all_cells().clear();
	    
	    finished_tracks.push_back(track1);
	    result_tracks.push_back(track1);
	    track1->get_end_scells().push_back(curr_cell);
	    track1->get_all_cells().push_back(curr_cell);
	  }	
	}
	
	
	//test
	// for (int j = 0; j!= result_tracks.size();j++){
	//   std::cout << result_tracks.at(j)->get_end_scells().at(0)->Get_Center().x/units::cm << " " 
	// 	    << result_tracks.at(j)->get_end_scells().at(1)->Get_Center().x/units::cm << " " 
	// 	    << std::endl;
	// }
	

      }
    }
    
  }
  

  return result_tracks;
}


WCTrackSelection WCVertex::BreakTracks(){
  WCTrackSelection result_tracks;
  if (tracks.size() > 1){
    float length = 0;
    
    for (int i=0;i!=tracks.size();i++){
       auto it = find(tracks.at(i)->get_end_scells().begin(),
		     tracks.at(i)->get_end_scells().end(),
		     msc);
       if (it != tracks.at(i)->get_end_scells().end()){
	 float length1 = fabs(tracks.at(i)->get_end_scells().at(0)->Get_Center().x-
			      tracks.at(i)->get_end_scells().at(1)->Get_Center().x)/units::cm;
	 if (length1 > length) length = length1;
       }
    }

    for (int i = 0;i!=tracks.size();i++){
      MergeSpaceCell *cell1 = tracks.at(i)->get_end_scells().at(0);
      MergeSpaceCell *cell2 = tracks.at(i)->get_end_scells().at(1);
      float dis1 = cell1->Get_Center().x - msc->Get_Center().x;
      float dis2 = cell2->Get_Center().x - msc->Get_Center().x;
      auto it = find(tracks.at(i)->get_end_scells().begin(),
		     tracks.at(i)->get_end_scells().end(),
		     msc);

      // std::cout << dis1 << " " << dis2 << " " << length << std::endl;

      if (it == tracks.at(i)->get_end_scells().end()){
	
	if (fabs(dis1)/units::cm > 0.65 &&
	    fabs(dis2)/units::cm > 0.65 && length > 2.0){
	  
	  WCTrack* primary_track = tracks.at(i);
	  WCTrack* secondary_track = new WCTrack(primary_track->get_mct());
	  if (cell1 == primary_track->replace_end_scells(msc)){
	    secondary_track->get_end_scells().clear();
	    secondary_track->get_end_scells().push_back(msc);
	    secondary_track->get_end_scells().push_back(cell1);
	  }else{
	    secondary_track->get_end_scells().clear();
	    secondary_track->get_end_scells().push_back(msc);
	    secondary_track->get_end_scells().push_back(cell2);
	  }
	  

	  tracks.push_back(secondary_track);

	  primary_track->ModifyCells();
	  secondary_track->ModifyCells();

	  result_tracks.push_back(primary_track);
	  result_tracks.push_back(secondary_track);

	  std::cout << "Break!" << " " << msc->Get_Center().x/units::cm << " " << 
	    msc->Get_Center().y/units::cm << " " << msc->Get_Center().z/units::cm << " " <<
	    fabs(cell1->Get_Center().x - msc->Get_Center().x)/units::cm << " " <<  fabs(cell2->Get_Center().x - msc->Get_Center().x)/units::cm << std::endl;
	  break;
	}

      }
    }
  }
  
  return result_tracks;
}

void WCVertex::ProcessTracks(WCTrackSelection& break_tracks){
  WCTrack* primary_track = break_tracks.at(0);
  WCTrack* secondary_track = break_tracks.at(1);

  // primary_track->ModifyCells();
  // secondary_track->ModifyCells();

  auto it = find(tracks.begin(),tracks.end(),primary_track);
  if (it != tracks.end()){
    float dis1 =  msc->Get_Center().x - primary_track->get_end_scells().at(0)->Get_Center().x;
    float dis2 =  msc->Get_Center().x - primary_track->get_end_scells().at(1)->Get_Center().x;
    
    if (dis1 * dis2 > 0){
      tracks.erase(it);
      tracks.push_back(secondary_track);
    }
    
  }
  
}


void WCVertex::OrganizeEnds(MergeSpaceCellSelection& cells1, int flag){
  WCTrackSelection removed;
  for (int i = 0;i!=tracks.size();i++){
    auto it = find(tracks.at(i)->get_end_scells().begin(),
  		   tracks.at(i)->get_end_scells().end(),
  		   msc);
    if (it == tracks.at(i)->get_end_scells().end()){
      MergeSpaceCell *cell1 = tracks.at(i)->get_end_scells().at(0);
      MergeSpaceCell *cell2 = tracks.at(i)->get_end_scells().at(1);
      if (fabs(cell1->Get_Center().x - msc->Get_Center().x)/units::cm > 0.65 &&
  	  fabs(cell2->Get_Center().x - msc->Get_Center().x)/units::cm > 0.65 && flag == 1){
  	removed.push_back(tracks.at(i));
      }else{
	if (flag == 1){
	  auto it1 = find(tracks.at(i)->get_all_cells().begin(),
			  tracks.at(i)->get_all_cells().end(),
			  msc);
	  float dis1 = tracks.at(i)->get_end_scells().at(0)->Get_Center().x - msc->Get_Center().x;
	  float dis2 = tracks.at(i)->get_end_scells().at(1)->Get_Center().x - msc->Get_Center().x;
	  if (it1 == tracks.at(i)->get_all_cells().end() && dis1*dis2 < 0){
	    removed.push_back(tracks.at(i));
	  }
	}
      }
    }
  }

  for (int i=0;i!=removed.size();i++){
    auto it = find(tracks.begin(),tracks.end(),removed.at(i));
    tracks.erase(it);
  }




  for (int i = 0; i!=tracks.size();i++){
    auto it = find(tracks.at(i)->get_end_scells().begin(),
  		   tracks.at(i)->get_end_scells().end(),
  		   msc);
    if (it == tracks.at(i)->get_end_scells().end()){
      tracks.at(i)->replace_end_scells(msc,&cells1);
      tracks.at(i)->ModifyCells();
    }
  }
}



bool WCVertex::CheckContain(MergeSpaceCell *cell){
  bool result = true;
  for (int i=0;i!=tracks.size();i++){
    WCTrack *track = tracks.at(i);
    //MergeClusterTrack& mct = track->get_mct();
    MergeSpaceCellSelection &cells = track->get_all_cells();
    auto it = find(cells.begin(),cells.end(),cell);
    if (it == cells.end()){
      result = false;
      break;
    }
  }
  return result;
}




int WCVertex::IsInside(WCVertex *vertex){
  int result = -1;
  
  // if (tracks.size()==1 && vertex->get_ntracks()==1){
  //   return -1;
  // }else if (tracks.size()>1 && vertex->get_ntracks()==1){
  //   return -1;
  // }else if (tracks.size()==1 && vertex->get_ntracks()>1){
  //   return -1;
  // }else if (tracks.size()>1 && vertex->get_ntracks()>1){
    result = 1;
    WCTrackSelection& temp_tracks = vertex->get_tracks();
    for (int i=0;i!=tracks.size();i++){
      WCTrack *track = tracks.at(i);
      auto it = find(temp_tracks.begin(),temp_tracks.end(),track);
      if (it == temp_tracks.end()){
	result = -1;
	break;
      }
    }
    
    if (result == 1){    
      for (int i = 0; i!=tracks.size();i++){
	WCTrack *track = tracks.at(i);
	auto it = find(track->get_all_cells().begin(),track->get_all_cells().end(),vertex->get_msc());
	if (it == track->get_all_cells().end()){
	  result = -1;
	  break;
	}
      }
    }
    
    
    if (result == 1 && tracks.size() == temp_tracks.size()){
      result = 0;
    }
  // }

  return result;
}

bool WCVertex::MergeVertex(WCVertex *vertex){
  if (msc != vertex->get_msc()){
    return false;
  }else{
    for (int i = 0; i!=vertex->get_tracks().size();i++){
      WCTrack *track = vertex->get_tracks().at(i);
      auto it = find(tracks.begin(),tracks.end(),track);
      if (it == tracks.end()){
	tracks.push_back(track);
	tracks_ky.push_back(vertex->get_ky(track));
	tracks_kz.push_back(vertex->get_kz(track));
      }
    }
    return true;
  }
}


bool WCVertex::AddVertex(WCVertex *vertex, int flag){
  bool result = false;
  
  // std::cout << "Xin1 " << std::endl;

  MergeSpaceCell *msc1 = msc;
  MergeSpaceCell *msc2 = vertex->get_msc();
  
  float dis = 7;
  


  if (fabs(msc1->Get_Center().x/units::mm-msc2->Get_Center().x/units::mm)<dis){
    // std::cout <<  fabs(msc1->Get_Center().x/units::mm-msc2->Get_Center().x/units::mm) << " " << msc1->Overlap(*msc2) << std::endl;
    if (msc1->Overlap(*msc2)){
      for (int i=0;i!=tracks.size();i++){
	WCTrack *track1 = tracks.at(i);
	auto it = find(vertex->get_tracks().begin(),vertex->get_tracks().end(),track1);
	if (it != vertex->get_tracks().end()){
	  result = true;
	  
	  //these two vertex share track1 ... 
	  if (flag == 1){
	    for (int j=0;j!=vertex->get_tracks().size();j++){
	      WCTrack *track = vertex->get_tracks().at(j);
	      auto it = find(track->get_all_cells().begin(),track->get_all_cells().end(),msc);
	      if (it == track->get_all_cells().end()){
		result = false;//track->Grow(msc);
	    	if (!result)
	    	  break;
	      }
	    }
	  } else if (flag==0){
	    //center must be at the end of common track
	    auto it1 = find(track1->get_end_scells().begin(),
	    		    track1->get_end_scells().end(),
	    		    msc);
	    if (it1 == track1->get_end_scells().end())
	      result = false;

	  }else if (flag == 2){
	    // No cut
	  }
	  

	  if (result){
	    // std::cout << tracks.size() << " " << tracks.at(0) << std::endl;
	    for (int j=0;j!=vertex->get_tracks().size();j++){
	      WCTrack *track2 = vertex->get_tracks().at(j);
	      //std::cout << track2 << std::endl;
	      auto it1 = find(tracks.begin(),tracks.end(),track2);
	      if (it1 == tracks.end()){
		tracks.push_back(track2);
	      }
	    }
	    // std::cout << tracks.size() << std::endl;
	    break;
	  }
	}
	
      }
    }
  }

  // std::cout << "Xin2 " << std::endl;

  return result;
}
