#ifndef WIRECELL_WIRECELLMAP
#define WIRECELL_WIRECELLMAP

#include "WireCellData/Cell.h"
#include "WireCellData/Wire.h"

#include <map>
#include <vector>

namespace WireCell {

    /// A data structure to associate a cell with some number of wires
    typedef std::map<const WireCell::Cell*, WireSelection> WireCellMap;
}

#endif
