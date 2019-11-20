#ifndef WIRECELL_WIRECHARGE_H
#define WIRECELL_WIRECHARGE_H

#include "WCPData/GeomWire.h"

namespace WCP {

    /// The Wire namespace holds types and functions related to a
    /// single charge on a wire.  The contents of this namespace have
    /// some symmetry with those of the Cell namespace.

    namespace Wire {

	/// An association between a geometry and GeomWire and a charge
	typedef std::pair<const GeomWire*,float> Charge;

	/// A simple collection of Wire::Charge objects.
	typedef std::vector<Charge> Group;

	/// A collection of Wire::Group objects
	typedef std::vector<Group> GroupCollection;

	/// Sum of charges
	double charge(const Group& group);

	/// Take a Wire::Group and return a Wire::GroupCollection which
	/// contains one Wire::Charge per group.
	GroupCollection singlets(const Group& group);

    } // namespace WCP::Wire
}     // namespace WCP


#endif
