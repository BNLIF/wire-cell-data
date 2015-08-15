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
    for (int j=0;j!=cells.size();j++){
      MergeSpaceCell *mscell = cells.at(j);
      double xc = mscell->Get_Center().x;
      
      int flag = 0;
      if (ntracks==1){
	if (fabs(xc-x0)/units::cm > 0.1 && fabs(xc-x0)/units::cm < 5)
	  flag = 1;
      }else if (ntracks>1){
	if (fabs(xc-x0)/units::cm > 0.7 && fabs(xc-x0)/units::cm < 5)
	  flag = 1;
	
      }
      
      if (flag == 1 ){
	// 1.0-5 cm
	
	//std::cout << xc/units::cm << " " << x0/units::cm << std::endl;

	SpaceCellSelection scells = mscell->Get_all_spacecell();
	for (int k=0;k!=scells.size();k++){
	  double x1 = scells.at(k)->x();
	  double y1 = scells.at(k)->y();
	  double z1 = scells.at(k)->z();
	  double q = scells.at(k)->q();

	  double a,b,c,d;
	  
	  a = y - ky[i]*x;
	  b = ky[i];
	  c = z - kz[i]*x;
	  d = kz[i];

	  double x2 = (-a*b-c*d+b*y1+d*z1+x1)/(1+b*b+d*d);
	  double y2 = a + b*x2;
	  double z2 = c + d*x2;

	  // std::cout << x1/units::cm << " " << x2/units::cm << " " << 
	  //   y1/units::cm << " " << y2/units::cm << " " << 
	  //   z1/units::cm << " " << z2/units::cm << " " <<  q << 
	  //   std::endl;

	  

	  // TVector3 n(1,ky[i],kz[i]);
	  // TVector3 p(x-x1,y-x1,z-x1);
	  // TVector3 v = p.Cross(n);
	  // double dis = v.Mag()/n.Mag();
	  
	  chi2 += (pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2))*q/0.15/0.15*3./units::cm/units::cm;
	  charge += q;
	}
      }
    }
  }

  //std::cout << chi2 << " " << charge << std::endl;
  
  return chi2/charge;
}

double MyFCN::operator() (const std::vector<double> & xx) const{
  return get_chi2(xx);
  
} 


int WCVertex::FindVertex(){



  MyFCN fcn(this);
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
  
  // // variable[0] = 62.0792 * units::cm;
  // // variable[1] = 35 * units::cm;
  // // variable[2] = 636*units::cm;

  // std::vector<double> variable1;
  // variable1.push_back(state.Value(0));
  // variable1.push_back(state.Value(1));
  // variable1.push_back(state.Value(2));
  // variable1.push_back(state.Value(3));

  // output
  //  std::cout << fcn.get_chi2(variable) << " Minimum: " << min.Fval() << " " << state.Value(0)/units::cm << " " <<state.Value(1)/units::cm << " " << state.Value(2)/units::cm << std::endl;
  // std::cout << "abc1: " << min.Fval() << " " << fcn.get_chi2(variable) << " " << fcn.get_chi2(variable1) << " " << state.Value(0)/units::cm << " " <<state.Value(1)/units::cm  << " " << state.Value(2) << " " << state.Value(3) << std::endl;

  //std::cout << fcn.get_chi2(variable) << std::endl;


  

 
   return 0;
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
    
    
  }
  


  
 
  
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
  	  fabs(cell2->Get_Center().x - msc->Get_Center().x)/units::cm > 0.65 ){
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
	    for (int j=0;j!=vertex->get_tracks().size();j++){
	      WCTrack *track2 = vertex->get_tracks().at(j);
	      auto it1 = find(tracks.begin(),tracks.end(),track2);
	      if (it1 == tracks.end()){
		tracks.push_back(track2);
	      }
	    }
	    break;
	  }
	}
	
      }
    }
  }

  // std::cout << "Xin2 " << std::endl;

  return result;
}
