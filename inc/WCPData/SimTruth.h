#ifndef WIRECELLDATA_SIMTRUTH
#define  WIRECELLDATA_SIMTRUTH

#include "WCPData/Point.h"

#include <set>

namespace WCP {

    /** WCP::SimTruth - information about simulation truth.
     */
    class SimTruth {
    public:
	SimTruth(float x, float y, float z, float q, int tdc, int trackid=-1);
	~SimTruth();		// do not inherit

	/// The X,Y,Z position of the deposited charge
	Point pos() const { return Point(_x,_y,_z); };
	/// The X position of the deposited charge
	float x() const { return _x; };
	/// The Y position of the deposited charge
	float y() const { return _y; };
	/// The Z position of the deposited charge
	float z() const { return _z; };
	/// The deposited charge
	float charge() const { return _q; };

	/// The TDC time bin where the majority of the deposited charge is collected.
	int tdc() const { return _tdc; };

	/// The simulation-specific ID of the track that produced the charge
	int trackid() const { return _trackid; };

    private:			// do not inherit
	float _x, _y, _z, _q;
	int _tdc, _trackid;

        friend std::ostream & operator<<(std::ostream &os, const SimTruth& st);
    };

    std::ostream & operator<<(std::ostream &os, const SimTruth& st);

    /// Compare SimTruth objects.
    struct SimTruthCompare {
	bool operator() (const SimTruth& a, const SimTruth& b) const {
	    if (a.tdc() < b.tdc()) {
		return true;
	    }
	    if (a.x() > b.x()) {
		return true;
	    }
	    return a.trackid() < b.trackid();
	}
    };


    /// An ordered, owning set of SimTruth objects.
    typedef std::set<WCP::SimTruth, SimTruthCompare> SimTruthSet;

    /// A selection on SimTruthSet.
    typedef std::vector<const SimTruth*> SimTruthSelection;
}

#endif
