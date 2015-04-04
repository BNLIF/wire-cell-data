#include "WireCellData/Slice.h"

using namespace WireCell;

Slice::Slice(int tbin, const Wire::Group& group)
    : _tbin(tbin)
    , _group(group)
{
}

Slice::~Slice()
{
}


void Slice::clear()
{
    _tbin = -1;
    _group.clear();
}

void Slice::reset(int tbin, const Wire::Group& group)
{
    _tbin = tbin;
    _group = group;
}
