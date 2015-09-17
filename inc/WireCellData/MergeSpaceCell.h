#ifndef MergeSpaceWireCellData_Cell_h
#define MergeSpaceWireCellData_Cell_h
#include "WireCellData/Point.h"
#include "WireCellData/SpaceCell.h"
#include "MergeGeomCell.h"
//#include "WireCellData/MergeClusterTrack.h"

#include <vector>
#include <map>
#include <list>
#include <set>

namespace WireCell {

    /** WireCell::Cell - information about one space cell
     */
    class MergeSpaceCell {
    public:
      MergeSpaceCell(){center_flag = 0;mcell = 0;};
      ~MergeSpaceCell();
	
      void AddSpaceCell(SpaceCell* cell){all_spacecell.push_back(cell);};
      SpaceCellSelection& Get_all_spacecell(){return all_spacecell;};
      
      bool Overlap(MergeSpaceCell& mcell, float num = 0.1);

      Point& Get_Center();
      float Get_Charge();

      float thickness(){return all_spacecell.front()->thickness();};

      bool CrossCell(Point &p, float theta, float phi, int flag = 0);
      double ClosestDis(Point &p);

      const MergeGeomCell* get_mcell(){return mcell;}
      void set_mcell(const MergeGeomCell *cell){mcell = cell;};

      void CalMinMax();
      double get_dy(){return fabs(max_y-min_y)/2.;};
      double get_dz(){return fabs(max_z-min_z)/2.;};
      double get_maxy(){return max_y;};
      double get_maxz(){return max_z;};
      double get_miny(){return min_y;};
      double get_minz(){return min_z;};


    protected:
      const MergeGeomCell* mcell;
      int center_flag;
      Point center;
      SpaceCellSelection all_spacecell;

      double max_y, min_y;
      double max_z, min_z;
      
};
    
    struct MergeSpaceCellCompare{
      bool operator() (MergeSpaceCell *a, MergeSpaceCell *b) const {
	
	if (a->Get_Center().x == b->Get_Center().x){
	  if (a->Get_all_spacecell().size()==b->Get_all_spacecell().size()){
	    return a > b;
	  }
	  return a->Get_all_spacecell().size() > b->Get_all_spacecell().size();
	}
	return a->Get_Center().x < b->Get_Center().x;
      }
    };

     

    /// Used to temporarily collect some subset
    typedef std::vector<WireCell::MergeSpaceCell*> MergeSpaceCellSelection;
    typedef std::map<MergeSpaceCell*, MergeSpaceCellSelection> MergeSpaceCellMap;
    typedef std::map<MergeSpaceCell*, MergeSpaceCell*> MergeSpaceCellMap1;
    typedef std::list<WireCell::MergeSpaceCell*> MergeSpaceCellList;
    typedef std::map<MergeSpaceCell*, int> MergeSpaceCellCounter;
    typedef std::set<MergeSpaceCell*, MergeSpaceCellCompare> MergeSpaceCellSet;

       
    // typedef std::map<MergeSpaceCell*, MergeClusterTrackSelection> MSC_MCT_Map;
    
}
#endif
