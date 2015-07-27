#include "WireCellData/MergeSpaceCell.h"

using namespace WireCell;

bool MergeSpaceCell::Overlap(MergeSpaceCell& mcell,float num){
  for (int i=0;i!=all_spacecell.size();i++){
    SpaceCell *cell1 = all_spacecell.at(i);
    for (int j=0;j!=mcell.Get_all_spacecell().size();j++){
      SpaceCell *cell2 = mcell.Get_all_spacecell().at(j);

      
      if (fabs(cell1->y()-cell2->y()) > 2.5*units::cm) break;
      if (fabs(cell1->z()-cell2->z()) > 2.5*units::cm) break;
      
      for (int i1=0;i1!=cell1->boundary().size();i1++){
	Point p = (cell1->boundary())[i1];
	for (int j1=0;j1!=cell2->boundary().size();j1++){
	  Point p1 = (cell2->boundary())[j1];
	  if (sqrt(pow(p.y-p1.y,2)+pow(p.z-p1.z,2))/units::m<0.003*num){
	    // std::cout << p.y << " " << p.z << " " << p1.y << " " << p1.z << std::endl;
	    return true;
	  }
	}
      }


    }
  }
  return false;
}

Point& MergeSpaceCell::Get_Center(){
  if (center_flag==1){
    return center;
  }else{
    center.x = 0;
    center.y = 0;
    center.z = 0;
    
    for (int i=0;i!=all_spacecell.size();i++){
      SpaceCell* cell = all_spacecell.at(i);
      center.x += cell->x();
      center.y += cell->y();
      center.z += cell->z();
    }
    
    center.x /= all_spacecell.size();
    center.y /= all_spacecell.size();
    center.z /= all_spacecell.size();
    
    return center;
  }
}
