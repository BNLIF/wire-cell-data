#ifndef WIRECELLDATA_TRACE_H
#define WIRECELLDATA_TRACE_H

#include <vector>
#include "TH1F.h"

namespace WCP {

    /// A disembodied ordered sequence of charge.
    typedef std::vector<float> ChargeSequence;
  
    /** Trace - a charge/time waveform signal on a wire.
      
      Note: zero-suppressed signals may be turned into an ordered
	collection of traces to save memory.

	See also Frame
     */

    struct Trace {
	/// The ID of the electronics channel from which this trace
	/// was read.  Note: use GeomDataSource to resolve wire id.
	int chid;
	
	/// The time bin relative to some absolute time at which the first ADC/charge exists
	int tbin;

	/// The charge in each time bin
	ChargeSequence charge;

        TH1F *hCharge; // same as charge except in the form of a histogram

        double channel_length;
      
    };

    typedef std::vector<Trace> TraceCollection;
}

#endif
