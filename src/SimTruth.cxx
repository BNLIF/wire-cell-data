#include "WCPData/SimTruth.h"

using namespace WCP;

SimTruth::SimTruth(float x, float y, float z, float q, int tdc, int trackid)
    : _x(x), _y(y), _z(z), _q(q), _tdc(tdc), _trackid(trackid)
{
}
SimTruth::~SimTruth()
{
}


