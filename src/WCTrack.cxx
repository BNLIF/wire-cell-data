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


MergeSpaceCell* WCTrack::replace_end_scells(MergeSpaceCell *cell2){
  Point p0 = cell2->Get_Center();
  MergeSpaceCell *cella = end_scells.at(0);
  MergeSpaceCell *cellb = end_scells.at(1);
  
  Point p1 = cella->Get_Center();
  Point p2 = cellb->Get_Center();
      
  float dis1 = sqrt(pow(p1.x-p0.x,2)+pow(p1.y-p0.y,2)+pow(p1.z-p0.z,2));
  float dis2 = sqrt(pow(p2.x-p0.x,2)+pow(p2.y-p0.y,2)+pow(p2.z-p0.z,2));

  
  std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm
	    << " " << p2.x/units::cm << " " << p2.y/units::cm <<  " " << p2.z/units::cm << " " << p0.x/units::cm << " " << p0.y/units::cm << " " << p0.z/units::cm << " " << dis1 << " " << dis2 << std::endl;

  if (dis1 > dis2){
    end_scells.pop_back();
    end_scells.push_back(cell2);
    return cellb;
  }else{
    MergeSpaceCell *cell1 = end_scells.at(1);
    end_scells.clear();
    end_scells.push_back(cell2);
    end_scells.push_back(cell1);
    return cella;
  }
  
  
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
