#include "WCPData/Frame.h"

using namespace WCP;

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
