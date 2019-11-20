#ifndef WIRECELLDATA_SLICE_H
#define WIRECELLDATA_SLICE_H

#include "WCPData/ChannelCharge.h"

#include <map>
#include <vector>


namespace WCP {

    /// A group of charges measured in channels at a given time bin
    class Slice {
    public:
	Slice(int tbin=-1, const Channel::Group& group = Channel::Group());
	~Slice();

	/// Forget the contents of the slice.
	void clear();

	/// Reset the values
	void reset(int tbin, const Channel::Group& group);
	
        /// Access the Channel::Group
        const Channel::Group& group() const { return _group; }
	Channel::Group& group() { return _group; }

	/// Access the associated time bin
	int tbin() const { return _tbin; }

    private:

	int _tbin;
	Channel::Group _group;

    };

}


#endif
