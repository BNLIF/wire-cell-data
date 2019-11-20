#ifdef __ROOTCLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link off all namespaces;
#pragma link C++ nestedclass;

#pragma link C++ namespace WCP;

#pragma link C++ namespace units;
#pragma link C++ defined_in namespace units;
#pragma link C++ global units::*;


#pragma link C++ class WCP::GeomCell;
#pragma link C++ class WCP::MergeGeomCell;
#pragma link C++ class WCP::MergeGeomCellSet;
#pragma link C++ class WCP::GeomCellMap;
#pragma link C++ class WCP::GeomCellSet;
#pragma link C++ class WCP::CellChargeMap;
#pragma link C++ class WCP::CellIndexMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WCP::GeomCellSet;
#pragma link C++ class WCP::GeomCellSelection;

#pragma link C++ class WCP::GeomWire;
#pragma link C++ class WCP::GeomWireSet;
#pragma link C++ class WCP::WireChargeMap;
#pragma link C++ class WCP::WireIndexMap;
#pragma link C++ class WCP::MergeGeomWire;
#pragma link C++ class WCP::GeomWireMap;
#pragma link C++ class WCP::GeomWireWireMap;
#pragma link C++ class WCP::GeomWireWiresMap;
// apparently, rootcling doesn't like sets
// #pragma link C++ class WCP::GeomWireSet;
#pragma link C++ class WCP::GeomWireSelection;

#pragma link C++ class WCP::MergeCellCluster;
#pragma link C++ class WCP::GeomCluster;
#pragma link C++ class WCP::GeomClusterSet;


#pragma link C++ class WCP::ChargeSequence;
#pragma link C++ class WCP::Trace;
#pragma link C++ class WCP::Frame;
//#pragma link C++ class WCP::WireCharge;
//#pragma link C++ class WCP::WireChargeCollection;
#pragma link C++ class WCP::Slice;

#pragma link C++ class WCP::Point;
#pragma link C++ class WCP::PointVector;

#pragma link C++ class WCP::Vector;
#pragma link C++ class WCP::VectorPair;

#pragma link C++ function WCP::box_intersection;

#endif
