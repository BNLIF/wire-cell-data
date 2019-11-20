#ifndef WIRECELL_WIRECELLMAP
#define WIRECELL_WIRECELLMAP

#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"

#include <map>

namespace WCP {

    /**
       Taken together, a (geometry) wire and a cell map provides a
       "binodal mesh" consisting of nodes of type GeomCell and
       Geomwire which connect to each other but not themselves.  Nodes
       hold pointers to the actual objects maintained elsewhere (see
       GeomCellSet and GeomWireSet).
     */

    typedef std::map<const GeomCell*, GeomWireSelection> GeomCellMap;
    typedef std::map<const GeomWire*, GeomCellSelection> GeomWireMap;
    typedef std::map<const GeomWire*, GeomCellList> GeomWireLMap;

    typedef std::map<const GeomWire*, GeomWire*> GeomWireWireMap;
    typedef std::map<const GeomWire*, GeomWireSelection> GeomWireWiresMap;

    typedef std::map<const GeomCell*, GeomCell*> GeomCellCellMap;
    typedef std::map<const GeomCell*, GeomCellSelection> GeomCellCellsMap;
}

#endif
