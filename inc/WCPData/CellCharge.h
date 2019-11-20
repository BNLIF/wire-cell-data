#ifndef WIRECELL_CELLCHARGE_H
#define WIRECELL_CELLCHARGE_H

#include "WCPData/GeomCell.h"
#include "WCPData/Point.h"

namespace WCP {

    /// The Cell namespace holds types and functions related to a
    /// single charge on a cell.  The contents of this namespace have
    /// some symmetry with those of the Wire namespace.

    namespace Cell {

	/// An association between a geometry and GeomCell and a charge
	typedef std::pair<const GeomCell*,float> Charge;

	/// A simple collection of Cell::Charge objects.
	typedef std::vector<Charge> Group;

	/// A collection of Cell::Group objects
	typedef std::vector<Group> GroupCollection;

	/// Sum of charges
	double charge(const Group& group);

	/// The charge weighted center of centers.
	WCP::Point center_of_charge(const Group& group);

	/// The total cross section
	double cross_section(const Group& group);

	/// Take a Cell::Group and return a Cell::GroupCollection
	/// which contains one Cell::Charge per group.
	GroupCollection singlets(const Group& group);

    } // namespace WCP::Cell
}     // namespace WCP


#endif
