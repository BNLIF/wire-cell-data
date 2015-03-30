#include "WireCellData/Frame.h"

using namespace WireCell;

Frame::Frame(int index, const TraceCollection& traces)
    : index(index)
    , traces(traces)
{
}
void Frame::clear()
{
    index = -1;
    traces.clear();
}
