#include "WireCellData/WireCharge.h"

using namespace WireCell;

double Wire::charge(const Wire::Group& group)
{
    double tot = 0;
    size_t nwires = group.size();
    for (size_t ind=0; ind<nwires; ++ind) {
	tot += group[ind].second;
    }
    return tot;
}


Wire::GroupCollection Wire::singlets(const Wire::Group& group)
{
    Wire::GroupCollection ret;
    size_t count = group.size();
    for (size_t ind=0; ind<count; ++ind) {
	Wire::Group g;
	g.push_back(group[ind]);
	ret.push_back(g);
    }
    return ret;
}

