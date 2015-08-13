#include "WireCellData/WCTrack.h"

using namespace WireCell;

WCTrack::WCTrack(MergeClusterTrack& mct)
  : mct(mct)
{
  MergeSpaceCellSelection& mcells = mct.Get_allmcells();
  end_scells.push_back(mcells.front());
  
  MergeSpaceCell* last_cell = mct.Get_LastMSCell();

  for (int i = mct.Get_allmcells().size()-1;i>=0; i--){
    MergeSpaceCell *cell = mct.Get_allmcells().at(i);
    if (cell->Get_Center().x!=last_cell->Get_Center().x){
      end_scells.push_back(mct.Get_allmcells().at(i+1));
      break;
    }
  }

}

WCTrack::~WCTrack(){
}


int WCTrack::TrackType(MergeSpaceCell& cell){
  int type = 0;
  // type == 1: short tracks
  // type == 2: straight tracks
  // type == 3: wiggle tracks 
  
  int time_length = mct.Get_TimeLength();
  if (time_length < 5){
    type = 1;
  }else{
    Point p = cell.Get_Center();
    
    //std::cout << theta << " " << phi << std::endl;
    int flag;
    Point p1 = mct.Get_FirstMSCell()->Get_Center();
    Point p2 = mct.Get_LastMSCell()->Get_Center();

    float dis1 = sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2));
    float dis2 = sqrt(pow(p.y-p2.y,2)+pow(p.z-p2.z,2));
    
    if (dis1 < dis2){
      flag = 1;
    }else{
      flag = -1;
    }
    

    for (int k=0;k!=5;k++){

      Point p3;
      if (flag==1){
	p3 =  mct.Get_MSCS(k).at(0)->Get_Center();
      }else{
	p3 =  mct.Get_MSCS(time_length-1-k).at(0)->Get_Center();
      }
    
      mct.SC_Hough(p3,p,10*units::cm,3);
      float theta = mct.Get_Theta();
      float phi = mct.Get_Phi();

      type = 2;
      
      for (int i=0;i!=5;i++){
	MergeSpaceCellSelection cells;
	if (flag == 1){
	  cells = mct.Get_MSCS(i);
	}else{
	  cells = mct.Get_MSCS(time_length-1-i);
	}
	int flag1 = 0;
	for (int j=0;j!=cells.size();j++){
	  MergeSpaceCell *cell = cells.at(j);
	  if (cell->CrossCell(p3,theta,phi)){
	    flag1 = 1;
	    break;
	  }
	}
	
	if (flag1==0){
	  type = 3;
	  break;
	}
	
      }
      
      if (type==2) break;
    }
  }
  
  return type;
}
