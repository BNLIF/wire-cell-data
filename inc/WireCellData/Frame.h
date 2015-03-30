#ifndef WIRECELL_FRAME_h
#define WIRECELL_FRAME_h

#include "WireCellData/Trace.h"

namespace WireCell {

    /** Frame - a collection of Trace objects read out over some
     * contemporaneous time period.
     */
    struct Frame {

	Frame(int index=-1, const TraceCollection& traces = TraceCollection());
	void clear();

	/// Index into collection of frames.
	int index;
	/// A collection of traces 
	TraceCollection traces;
    };

};	

#endif
