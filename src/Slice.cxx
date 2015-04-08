#include "WireCellData/Slice.h"

using namespace WireCell;

Slice::Slice(int tbin, const Channel::Group& group)
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

void Slice::reset(int tbin, const Channel::Group& group)
{
    _tbin = tbin;
    _group = group;
}
