#ifndef MergeSpaceWireCellData_Cell_h
#define MergeSpaceWireCellData_Cell_h

#include "WireCellData/SpaceCell.h"

#include <set>
#include <vector>
#include <iostream>
#include <list>

namespace WireCell {

    /** WireCell::Cell - information about one space cell
     */
    class MergeSpaceCell {
    public:
      MergeSpaceCell(){};
      ~MergeSpaceCell(){};
	
      void AddSpaceCell(SpaceCell* cell){all_spacecell.push_back(cell);};
      SpaceCellSelection& Get_all_spacecell(){return all_spacecell;};
	
    protected:
      
      SpaceCellSelection all_spacecell;
      
};
    
    /// Used to temporarily collect some subset
    typedef std::vector<WireCell::MergeSpaceCell*> MergeSpaceCellSelection;
}
#endif
