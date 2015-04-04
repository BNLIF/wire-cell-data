#ifndef WIRECELLDATA_SLICE_H
#define WIRECELLDATA_SLICE_H

#include "WireCellData/WireCharge.h"

#include <map>
#include <vector>


namespace WireCell {

    /// A group of charges on wires at a given time bin
    class Slice {
    public:
	Slice(int tbin=-1, const Wire::Group& group = Wire::Group());
	~Slice();

	/// Forget the contents of the slice.
	void clear();

	/// Reset the values
	void reset(int tbin, const Wire::Group& group);
	
        /// Access the Wire::Group
        const Wire::Group& group() const { return _group; }
	Wire::Group& group() { return _group; }

	/// Access the associated time bin
	int tbin() const { return _tbin; }

    private:

	int _tbin;
	Wire::Group _group;

    };

}


#endif
