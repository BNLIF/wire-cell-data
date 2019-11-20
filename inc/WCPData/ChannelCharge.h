#ifndef WIRECELL_CHANNELCHARGE_H
#define WIRECELL_CHANNELCHARGE_H

#include <map>			// for pair
#include <vector>

namespace WCP {

    /// The Channel namespace holds types and functions related to a
    /// single charge on a channel.  The contents of this namespace
    /// have some symmetry with those of the Cell namespace.

    namespace Channel {

	/// An association between a channel and a charge
	typedef std::pair<int, float> Charge;

	/// A simple collection of Channel::Charge objects.
	typedef std::vector<Charge> Group;

	/// A collection of Channel::Group objects
	typedef std::vector<Group> GroupCollection;

	/// Sum of charges
	double charge(const Group& group);

	/// Take a Channel::Group and return a
	/// Channel::GroupCollection which contains one
	/// Channel::Charge per group.
	GroupCollection singlets(const Group& group);

    }
}    


#endif
