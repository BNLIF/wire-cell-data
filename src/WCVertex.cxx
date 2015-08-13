#include "WireCellData/WCVertex.h"

using namespace WireCell;

WCVertex::WCVertex(MergeSpaceCell& msc)
  : msc(msc)
{
  center = msc.Get_Center();
}

WCVertex::~WCVertex(){
}


Point WCVertex::Center(){
  
  return center;
}

void WCVertex::Add(WCTrack* track){
  tracks.push_back(track);
}

void WCVertex::OrganizeTracks(){
  
  if (tracks.size()>1){
    
    // find the track where the vertex is one of the end track
    WCTrack * end_track;
    for (int i = 0;i!=tracks.size();i++){
      auto it = find(tracks.at(i)->get_end_scells().begin(),
		     tracks.at(i)->get_end_scells().end(),
		     &msc);
      if (it != tracks.at(i)->get_end_scells().end()){
	end_track = tracks.at(i);
	break;
      }
    }
    
    MergeClusterTrack& mct = end_track->get_mct();
    //find out where is the mct
    int flag;
    Point p0 = msc.Get_Center();
    Point p1 = mct.Get_FirstMSCell()->Get_Center();
    Point p2 = mct.Get_LastMSCell()->Get_Center();
    float dis1 = sqrt(pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2));
    float dis2 = sqrt(pow(p0.x-p2.x,2)+pow(p0.y-p2.y,2)+pow(p0.z-p2.z,2));
    
    if (dis1 < dis2){
      flag = 1;
    }else{
      flag = -1;
    }
    
    
    int time_length = mct.Get_TimeLength();
    MergeSpaceCell* nvertex = &msc;

    //std::cout << time_length << " " << flag << " " << nvertex->Get_Center().x/units::cm << std::endl;
    if (flag ==1){
      int flag_crawl = 0;      
      for (int i=0;i!=time_length;i++){
	MergeSpaceCellSelection& cells = mct.Get_MSCS(i);
	int flag_common = 0;
	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *tcell = cells.at(j);
	  if (tcell == &msc){
	    flag_crawl = 1; // start to crawl
	  }
	  if (flag_crawl == 1){
	    //std::cout << this->CheckContain(tcell) << std::endl;
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
	  }
	}
	if (flag_crawl == 1 && flag_common == 0)
	  break;
      }
    }else{
      int flag_crawl = 0;      
      for (int i=0;i!=time_length;i++){
	MergeSpaceCellSelection& cells = mct.Get_MSCS(time_length-1-i);
	int flag_common = 0;
	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *tcell = cells.at(j);
	  if (tcell == &msc){
	    flag_crawl = 1; // start to crawl
	  }
	  if (flag_crawl == 1){
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
	  }
	}
	if (flag_crawl == 1 && flag_common == 0)
	  break;
      }
    }
    msc = *nvertex;
    center = msc.Get_Center();
    //    std::cout << msc.Get_Center().x/units::cm << " " << nvertex->Get_Center().x/units::cm << std::endl;
    
  }
  

// tracks.at(0)->replace_end_scells(&msc);
//   }else{
  // for (int i = 0; i!=tracks.size();i++){
  //   auto it = find(tracks.at(i)->get_end_scells().begin(),
  // 		   tracks.at(i)->get_end_scells().end(),
  // 		   &msc);
  //   if (it == tracks.at(i)->get_end_scells().end()){
  //     //    tracks.at(i)->replace_end_scells(&msc);
  //   }
  // }
 
  
}

bool WCVertex::CheckContain(MergeSpaceCell *cell){
  bool result = true;
  for (int i=0;i!=tracks.size();i++){
    WCTrack *track = tracks.at(i);
    MergeClusterTrack& mct = track->get_mct();
    MergeSpaceCellSelection &cells = mct.Get_allmcells();
    auto it = find(cells.begin(),cells.end(),cell);
    if (it == cells.end()){
      result = false;
      break;
    }
  }
  return result;
}




int WCVertex::IsInside(WCVertex *vertex){
  int result = 1;
  WCTrackSelection& temp_tracks = vertex->get_tracks();
  for (int i=0;i!=tracks.size();i++){
    WCTrack *track = tracks.at(i);
    auto it = find(temp_tracks.begin(),temp_tracks.end(),track);
    if (it == temp_tracks.end()){
      result = -1;
      break;
    }
  }
  
  if (result == 1 && tracks.size() == temp_tracks.size()){
    result = 0;
  }
  

  return result;
}


bool WCVertex::AddVertex(WCVertex *vertex){
  bool result = false;
  
  MergeSpaceCell *msc1 = &msc;
  MergeSpaceCell *msc2 = vertex->get_msc();
  
  if (fabs(msc1->Get_Center().x/units::mm-msc2->Get_Center().x/units::mm)<5){
    // std::cout <<  fabs(msc1->Get_Center().x/units::mm-msc2->Get_Center().x/units::mm) << " " << msc1->Overlap(*msc2) << std::endl;
    if (msc1->Overlap(*msc2)){
      for (int i=0;i!=tracks.size();i++){
	WCTrack *track1 = tracks.at(i);
	auto it = find(vertex->get_tracks().begin(),vertex->get_tracks().end(),track1);
	if (it != vertex->get_tracks().end()){
	  result = true;
	  
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

  return result;
}
