#include "WCPData/GeomWire.h"

#include <cassert>
#include <iostream>

using namespace std;
using namespace WCP;

int main()
{
    // these do not imply actual or required numbering conventions.
    // ident, wire plane type, wire plane index, channel, 
    GeomWireSet gws;

    // construct out of order intentionally so check by plane+index
    gws.insert(GeomWire(111, kYwire, 3, 303));
    gws.insert(GeomWire(222, kVwire, 0, 300));
    gws.insert(GeomWire(333, kUwire, 1, 301));
    gws.insert(GeomWire(444, kYwire, 2, 302));

    GeomWireSet::iterator it, done = gws.end();

    for (it = gws.begin(); it != done; ++it) {
	cerr << *it << endl;
    }

    int *want, wanted[] = {333,222,444,111};
    int nerrors = 0;
    for (want=wanted, it = gws.begin(); it != done; ++it, ++want) {
	int got = it->ident();
	if (got != *want) {
	    cerr << "got:" << got << " want:" << *want << endl;
	    ++nerrors;
	}
    }
    assert (!nerrors);

    return 0;
}
