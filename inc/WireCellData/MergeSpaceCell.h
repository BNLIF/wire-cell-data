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
	
      void AddSpaceCell(SpaceCell* cell); 
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

      GeomWireSelection get_uwires() const { if(mcell==0){
	  return uwires;
	}else{
	  return mcell->get_uwires();
	}
      };
      GeomWireSelection get_vwires() const { if(mcell==0){
	  return vwires;
	}else{
	  return mcell->get_vwires();
	}
      };
      GeomWireSelection get_wwires() const {if (mcell==0){
	  return wwires;
	}else{
	  return mcell->get_wwires();
	}
      };

      void set_id(int num){mcell_id = num;};
      int get_id(){return mcell_id;};

    protected:
      int mcell_id;
      const MergeGeomCell* mcell;
      int center_flag;
      Point center;
      GeomWireSelection uwires;
      GeomWireSelection vwires;
      GeomWireSelection wwires;
      SpaceCellSelection all_spacecell;

      double max_y, min_y;
      double max_z, min_z;
      
};
    
    struct MergeSpaceCellCompare{
      bool operator() (MergeSpaceCell *a, MergeSpaceCell *b) const {
	
	if (fabs(a->Get_Center().x - b->Get_Center().x) < 0.1*units::mm ){
	  if (a->Get_all_spacecell().size()==b->Get_all_spacecell().size()){
	    return a > b;
	  }
	  return a->Get_all_spacecell().size() > b->Get_all_spacecell().size();
	}
	return a->Get_Center().x < b->Get_Center().x;
      }
    };

    struct MergeSpaceCellCompare1{
      bool operator() (MergeSpaceCell *a, MergeSpaceCell *b) const {
	
	if (fabs(a->Get_Center().x - b->Get_Center().x)<0.1*units::mm){
	  if (a->Get_all_spacecell().size()==b->Get_all_spacecell().size()){
	    return a > b;
	  }
	  return a->Get_all_spacecell().size() > b->Get_all_spacecell().size();
	}
	return a->Get_Center().x > b->Get_Center().x;
      }
    };

     

    /// Used to temporarily collect some subset
    typedef std::vector<WireCell::MergeSpaceCell*> MergeSpaceCellSelection;
    typedef std::map<MergeSpaceCell*, MergeSpaceCellSelection> MergeSpaceCellMap;
    typedef std::map<MergeSpaceCell*, MergeSpaceCell*> MergeSpaceCellMap1;
    typedef std::list<WireCell::MergeSpaceCell*> MergeSpaceCellList;
    typedef std::map<MergeSpaceCell*, int> MergeSpaceCellCounter;
    typedef std::set<MergeSpaceCell*, MergeSpaceCellCompare> MergeSpaceCellSet;
    typedef std::set<MergeSpaceCell*, MergeSpaceCellCompare1> MergeSpaceCellSet1;

       
    // typedef std::map<MergeSpaceCell*, MergeClusterTrackSelection> MSC_MCT_Map;
    
}
#endif
