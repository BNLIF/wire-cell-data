#include "WireCellData/MergeSpaceCell.h"

using namespace WireCell;

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
