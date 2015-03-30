#ifndef WIRECELLDATA_SLICE_H
#define WIRECELLDATA_SLICE_H

#include <map>
#include <vector>

namespace WireCell {

    /// An association of wire ID and charge
    typedef std::pair<int,float> WireCharge;

    /// A collection of wire charges
    typedef std::vector<WireCharge> WireChargeCollection;

    /// A collection of wire charges at a given time bin
    struct Slice {
	Slice(int tbin=-1, const WireChargeCollection& charge = WireChargeCollection());
	void clear();

	int tbin;
	WireChargeCollection charge;

    };

}


#endif
