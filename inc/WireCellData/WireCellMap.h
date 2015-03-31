#ifndef WIRECELL_WIRECELLMAP
#define WIRECELL_WIRECELLMAP

#include "WireCellData/Cell.h"
#include "WireCellData/Wire.h"

#include <map>

namespace WireCell {

    /**
       Taken together, a wire and a cell map is a mesh consisting of
       two types of nodes (cell and wire) which connect to each other
       but not themselves.  Nodes hold pointers to the actual objects
       maintained elsewhere (see CellSet and WireSet).
     */

    typedef std::map<const Cell*, WireSelection> CellMap;
    typedef std::map<const Wire*, CellSelection> WireMap;
    
}

#endif
