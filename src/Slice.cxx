#include "WireCellData/Slice.h"

using namespace WireCell;

Slice::Slice(int tbin, const WireChargeCollection& charge)
    : tbin(tbin)
    , charge(charge)
{
}

void Slice::clear()
{
    tbin = -1;
    charge.clear();
}
