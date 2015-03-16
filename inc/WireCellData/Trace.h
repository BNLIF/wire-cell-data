#ifndef WIRECELLDATA_TRACE_H
#define WIRECELLDATA_TRACE_H

#include <vector>

namespace WireCellData {

    /// A disembodied ordered sequence of charge.
    typedef std::vector<float> ChargeSequence;

    /**
       Hold a trace or signal as a function of time bins on one wire.

       A trace is a portion of possibly a larger signal (eg, a segment
       which is above some threshold or separated into frames).  
     */

    struct Trace {
	/// The ID number of the wire on which the trace was measured.
	int wid;
	
	/// The time bin relative to some absolute time at which the first ADC/charge exists
	int tbin;

	/// The charge in each time bin
	ChargeSequence charge;
    };

    /// A collection of traces, such as would be read out in a single
    /// frame.  Multiple Traces from the same wire may be present in
    /// the same frame (eg, due to zero suppression)
    typedef std::vector<Trace> Frame;
}

#endif
