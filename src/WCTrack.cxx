#include "WireCellData/WCTrack.h"
#include "WireCellData/Plane.h"

using namespace WireCell;

bool WCTrack::fine_tracking(Point &p1, double ky1, double kz1, Point &p2, double ky2, double kz2){
  //if (fine_tracking_flag==1) return false;
  fine_tracking_flag = 1;
  
  centerVP.clear();
  PointVector frontVP;
  PointVector backVP;

  //initial filling of the center vector
  for (int i=0; i!=all_cells.size();i++){
    Point p;
    p.x = all_cells.at(i)->Get_Center().x;
    p.y = all_cells.at(i)->Get_Center().y;
    p.z = all_cells.at(i)->Get_Center().z;
    
    //std::cout << p.x/units::cm << " " << p.y/units::cm << " " << p.z/units::cm << std::endl;
    
    if (i==0){
      centerVP.push_back(p);
      frontVP.push_back(p);
      backVP.push_back(p);
      centerVP_cells.push_back(all_cells.at(i));
    }else{
      if (fabs(p.x-centerVP.at(centerVP.size()-1).x)>0.1*units::mm){
     	centerVP.push_back(p);
    	frontVP.push_back(p);
    	backVP.push_back(p);
    	centerVP_cells.push_back(all_cells.at(i));
      }else{
    	if (all_cells.at(i)->Get_all_spacecell().size() > centerVP_cells.at(centerVP_cells.size()-1)->Get_all_spacecell().size()){
    	  centerVP.at(centerVP_cells.size()-1) = p;
    	  frontVP.at(centerVP_cells.size()-1) = p;
    	  backVP.at(centerVP_cells.size()-1) = p;
    	  centerVP_cells.at(centerVP_cells.size()-1) = all_cells.at(i);
    	}
      }
    }
  }

  // judge if first cell is closer to p1 or p2
  Point pf;
  pf.x = centerVP_cells.at(0)->Get_Center().x;
  pf.y = centerVP_cells.at(0)->Get_Center().y;
  pf.z = centerVP_cells.at(0)->Get_Center().z;
  float dis1 = pow(pf.x-p1.x,2);
  float dis2 = pow(pf.x-p2.x,2);
  int order;
  if (dis1<dis2){
    order = 1;  // firt cell closer to p1, 
  }else{
    order = 2;  // first cell closer to p2,
  }

  // iteration procedure ... 
  for (int k=0;k!=3;k++){
    
    // if (order == 1){
	  
    for (int i=0;i!=centerVP.size();i++){
      if (fabs(centerVP.at(i).x - p1.x) < 0.65*units::cm){ // twice the time difference
	double min_dis = 1e9;
	Point min_point;
	Point p3;
	p3.x = p1.x + 1;
	p3.y = p1.y + ky1;
	p3.z = p1.z + kz1;
	Line l1(p1,p3);

	for (int j=0;j!=centerVP_cells.at(i)->Get_all_spacecell().size();j++){
	  SpaceCell *cell = centerVP_cells.at(i)->Get_all_spacecell().at(j);
	  Point p;
	  p.x = cell->x();
	  p.y = cell->y();
	  p.z = cell->z();
	  double dis = l1.closest_dis(p);
	  if (dis < min_dis){
	    min_dis = dis;
	    min_point = p;
	  }
	}
	
	if (min_dis > 3*units::mm){
	  frontVP.at(i) = centerVP.at(i);
	  backVP.at(i) = centerVP.at(i);
	}else{
	  centerVP.at(i) = min_point;
	  frontVP.at(i) = centerVP.at(i);
	  backVP.at(i) = centerVP.at(i);
	}
	
      }else if (fabs(centerVP.at(i).x-p2.x) < 0.65*units::cm){ // twice the time difference
	double min_dis = 1e9;
	Point min_point;
	Point p3;
	p3.x = p2.x + 1;
	p3.y = p2.y + ky2;
	p3.z = p2.z + kz2;
	Line l1(p1,p3);

	for (int j=0;j!=centerVP_cells.at(i)->Get_all_spacecell().size();j++){
	  SpaceCell *cell = centerVP_cells.at(i)->Get_all_spacecell().at(j);
	  Point p;
	  p.x = cell->x();
	  p.y = cell->y();
	  p.z = cell->z();
	  double dis = l1.closest_dis(p);
	  if (dis < min_dis){
	    min_dis = dis;
	    min_point = p;
	  }
	}
	
	if (min_dis > 3*units::mm){
	  frontVP.at(i) = centerVP.at(i);
	  backVP.at(i) = centerVP.at(i);
	}else{
	  centerVP.at(i) = min_point;
	  frontVP.at(i) = centerVP.at(i);
	  backVP.at(i) = centerVP.at(i);
	}
      }else{
	// go through each of the merge blob
	Plane plane1(centerVP.at(i),p1,p2);
	if (plane1.sameline()){
	  plane1.get_p1().y += 0.1*units::mm; // move the center a little bit
	}
	//get the plane
	Point pp1=centerVP.at(i);
	Point pp2=centerVP.at(i);
	pp1.y = pp1.y+1*units::cm;
	pp2.z = pp2.z+1*units::cm;
	Plane plane2(centerVP.at(i),pp1,pp2);
	//find the intersection line
	Line l1(p1,p2);
	Line& l2 = plane1.CrossLineCommonPoint(plane2);
	TVector3& v1 = l1.vec();
	TVector3& v2 = l2.vec();
	if (v1.Dot(v2)<0){
	  l2.ReverseDir();
	}
	
	//Now loop through all the cells to find the one satisfy the cuts
	for (int j=0;j!=centerVP_cells.at(i)->Get_all_spacecell().size();j++){
	  SpaceCell *cell = centerVP_cells.at(i)->Get_all_spacecell().at(j);
	  Point p;
	  p.x = cell->x();
	  p.y = cell->y();
	  p.z = cell->z();
	  double dis = l2.closest_dis(p);
	  if (dis < 4.5 * units::mm){
	    TVector3 v3(p.x-centerVP.at(i).x,p.y-centerVP.at(i).y,p.z-centerVP.at(i).z);
	    TVector3 v4(frontVP.at(i).x-centerVP.at(i).x,frontVP.at(i).y-centerVP.at(i).y,frontVP.at(i).z-centerVP.at(i).z);
	    TVector3 v5(backVP.at(i).x-centerVP.at(i).x,backVP.at(i).y-centerVP.at(i).y,backVP.at(i).z-centerVP.at(i).z);
	    
	    double dist1 = v1.Dot(v3);
	    double dist2 = v1.Dot(v4);
	    double dist3 = v1.Dot(v5);
	    
	    if (dist1 > dist2){
	      frontVP.at(i).x = p.x;
	      frontVP.at(i).y = p.y;
	      frontVP.at(i).z = p.z;
	    }else if (dist1 < dist3){
	      backVP.at(i).x = p.x;
	      backVP.at(i).y = p.y;
	      backVP.at(i).z = p.z;
	    }
	  }
	}
      }
    } // for loop
    
    //} // end order ... 
  
    // //Now recalculate the center
    
    for (int i=0;i!=centerVP.size();i++){
    // //   if (centerVP.size()==1){
    // // 	centerVP.at(i).x = (frontVP.at(i).x + backVP.at(i).x)/2.;
      centerVP.at(i).y = (frontVP.at(i).y + backVP.at(i).y)/2.;
      centerVP.at(i).z = (frontVP.at(i).z + backVP.at(i).z)/2.;
    // //   }else{
    // // 	if (i==0){
    // // 	}else if (i==centerVP.size()-1){
    // // 	}else{
    // // 	}
    // //   }
    }

  }

  // construct the differential angles ... 
  for (int i=0;i!=centerVP.size()-1;i++){
    TVector3 dir(centerVP.at(i+1).x-centerVP.at(i).x,
		 centerVP.at(i+1).y-centerVP.at(i).y,
		 centerVP.at(i+1).z-centerVP.at(i).z);
    centerVP_theta.push_back(dir.Theta());
    centerVP_phi.push_back(dir.Phi());
    if (i==centerVP.size()-2){
      centerVP_theta.push_back(dir.Theta());
      centerVP_phi.push_back(dir.Phi());
    }
  }
  


  return true;

}

WCTrack::WCTrack(MergeClusterTrack& mct)
  : mct(mct)
{
  fine_tracking_flag = 0;
  MergeSpaceCellSelection& mcells = mct.Get_allmcells();
  end_scells.push_back(mcells.front());
  
  MergeSpaceCell* last_cell = mct.Get_LastMSCell();
  // end_scells.push_back(last_cell);

  for (int i = mct.Get_allmcells().size()-1;i>=0; i--){
    MergeSpaceCell *cell = mct.Get_allmcells().at(i);
    if (cell->Get_Center().x!=last_cell->Get_Center().x){
      end_scells.push_back(mct.Get_allmcells().at(i+1));
      break;
    }
  }

  for (int i=0;i!=mcells.size();i++){
    all_cells.push_back(mcells.at(i));
  }


}

WCTrack::~WCTrack(){
}

bool WCTrack::Grow(MergeSpaceCell *cell, int flag){
  MergeSpaceCell *cell1 = end_scells.at(0);
  MergeSpaceCell *cell2 = end_scells.at(1);

  //test with the first cell
  if ((cell->Get_Center().x - cell1->Get_Center().x)/units::cm >0.4 && cell->Overlap(*cell1)){
    if (flag == 0){
      end_scells.at(0) = cell;
      all_cells.insert(all_cells.begin(),cell);
    }
    return true;
  }

  //test with the last cell
  if ((cell->Get_Center().x - cell2->Get_Center().x)/units::cm <0.4 && cell->Overlap(*cell2)){
    if (flag == 0){
      end_scells.at(1) = cell;
      all_cells.push_back(cell);    
    }
    return true;
  }

  return false;
}


void WCTrack::ModifyCells(){
  MergeSpaceCellSelection temp;
  MergeSpaceCell *cell1 = end_scells.at(0);
  MergeSpaceCell *cell2 = end_scells.at(1);

  for (int i=0;i!=all_cells.size();i++){
    MergeSpaceCell *cell = all_cells.at(i);
    if ((cell->Get_Center().x >= cell1->Get_Center().x && 
	 cell->Get_Center().x <= cell2->Get_Center().x) ||
	(cell->Get_Center().x <= cell1->Get_Center().x && 
	 cell->Get_Center().x >= cell2->Get_Center().x)){
      temp.push_back(cell);
    }
  }
  all_cells.clear();
  all_cells = temp;
  
}


MergeSpaceCell* WCTrack::replace_end_scells(MergeSpaceCell *cell2, MergeSpaceCellSelection* cells){
  Point p0 = cell2->Get_Center();
  MergeSpaceCell *cella = end_scells.at(0);
  MergeSpaceCell *cellb = end_scells.at(1);
  
  Point p1 = cella->Get_Center();
  Point p2 = cellb->Get_Center();
      
  float dis1 = sqrt(pow(p1.x-p0.x,2)+pow(p1.y-p0.y,2)+pow(p1.z-p0.z,2));
  float dis2 = sqrt(pow(p2.x-p0.x,2)+pow(p2.y-p0.y,2)+pow(p2.z-p0.z,2));
  
  float dis10 = fabs(p1.x-p0.x);
  float dis20 = fabs(p2.x-p0.x);
  
  // std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm
  // 	    << " " << p2.x/units::cm << " " << p2.y/units::cm <<  " " << p2.z/units::cm << " " << p0.x/units::cm << " " << p0.y/units::cm << " " << p0.z/units::cm << " " << dis1 << " " << dis2 << std::endl;



  if (dis10 > dis20){
    if (cells!=0){
      auto it = find(cells->begin(),cells->end(),cellb);
      if (it!=cells->end())
	return 0;
    }
    end_scells.pop_back();
    end_scells.push_back(cell2);
    return cellb;
  }else if (dis10 < dis20){
    if (cells!=0){
      auto it = find(cells->begin(),cells->end(),cella);
      if (it!=cells->end())
	return 0;
    }
    MergeSpaceCell *cell1 = end_scells.at(1);
    end_scells.clear();
    end_scells.push_back(cell2);
    end_scells.push_back(cell1);
    return cella;
  }else{
    if (dis1 > dis2){
      if (cells!=0){
	auto it = find(cells->begin(),cells->end(),cellb);
	if (it!=cells->end())
	  return 0;
      }
      end_scells.pop_back();
      end_scells.push_back(cell2);
      return cellb;
    }else{
      if (cells!=0){
	auto it = find(cells->begin(),cells->end(),cella);
	if (it!=cells->end())
	  return 0;
      }
      MergeSpaceCell *cell1 = end_scells.at(1);
      end_scells.clear();
      end_scells.push_back(cell2);
      end_scells.push_back(cell1);
      return cella;
    }
  }
  
  
}


void WCTrack::ReplaceEndCell(MergeSpaceCell *cell1, MergeSpaceCell *cell2){
  if (end_scells.at(0) == cell1){
    end_scells.at(0) = cell2;
  }
  if (end_scells.at(1) == cell1){
    end_scells.at(1) = cell2;
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
