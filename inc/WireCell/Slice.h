#ifndef WIRECELL_SLICE_H
#define WIRECELL_SLICE_H

#include <map>
#inlcude <vector>

namespace WireCell {

    /// An association of wire ID and charge
    typedef std::pair<int,float> WireCharge;

    /// A collection of wire charges
    typedef std::vector<WireCharge> WireChargeCollection;

    /// A collection of wire charges at a given time bin
    struct Slice {
	int tbin;
	WireChargeCollection charge;
    };

}


#endif
