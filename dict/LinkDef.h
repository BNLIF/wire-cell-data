#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace units;
#pragma link C++ namespace WireCell;

#pragma link C++ class WireCell::GeomCell;
#pragma link C++ class WireCell::GeomCellMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WireCell::GeomCellSet;
#pragma link C++ class WireCell::GeomCellSelection;

#pragma link C++ class WireCell::GeomWire;
#pragma link C++ class WireCell::GeomWireMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WireCell::GeomWireSet;
#pragma link C++ class WireCell::GeomWireSelection;

#pragma link C++ class WireCell::ChargeSequence;
#pragma link C++ class WireCell::Trace;
#pragma link C++ class WireCell::Frame;
#pragma link C++ class WireCell::WireCharge;
#pragma link C++ class WireCell::WireChargeCollection;
#pragma link C++ class WireCell::Slice;

#pragma link C++ class WireCell::Point;
#pragma link C++ class WireCell::PointVector;
#endif
