#ifndef MergeSpaceWireCellData_Cell_h
#define MergeSpaceWireCellData_Cell_h
#include "WireCellData/Point.h"
#include "WireCellData/SpaceCell.h"


#include <vector>
#include <map>

namespace WireCell {

    /** WireCell::Cell - information about one space cell
     */
    class MergeSpaceCell {
    public:
      MergeSpaceCell(){center_flag = 0;};
      ~MergeSpaceCell(){};
	
      void AddSpaceCell(SpaceCell* cell){all_spacecell.push_back(cell);};
      SpaceCellSelection& Get_all_spacecell(){return all_spacecell;};
      
      bool Overlap(MergeSpaceCell& mcell, float num = 0.1);

      Point& Get_Center();

      double thickness(){return all_spacecell.front()->thickness();};

    protected:
      int center_flag;
      Point center;
      SpaceCellSelection all_spacecell;
      
};
    
    /// Used to temporarily collect some subset
    typedef std::vector<WireCell::MergeSpaceCell*> MergeSpaceCellSelection;
    typedef std::map<MergeSpaceCell*, MergeSpaceCellSelection> MergeSpaceCellMap;
    typedef std::map<MergeSpaceCell*, int> MergeSpaceCellCounter;
    
    
}
#endif
