#include "WireCellData/MergeGeomWire.h"

#include <vector>
#include <cmath>
using namespace std;
using namespace WireCell;


MergeGeomWire::MergeGeomWire(int ident, GeomWireSelection wires){
  
  _ident = ident;
  _plane = wires.at(0)->plane();
  
  _index = -1;
  _channel = -1;
  _point1 = Point();
  _point2 = Point();
  
  time_slice = -1;
  
  for(int i=0;i!=wires.size();i++){
    wire_all.push_back(wires.at(i));
  }
}

MergeGeomWire::MergeGeomWire(int ident, const WireCell::GeomWire& wire)
{
  _ident = ident;
  _plane = wire.plane();
  
  _index = -1;
  _channel = -1;
  _point1 = Point();
  _point2 = Point();
  
  time_slice = -1;

  wire_all.push_back(&wire);
  
}

MergeGeomWire::MergeGeomWire(const WireCell::MergeGeomWire& wire)
{
  _ident = wire.ident();
  _plane = wire.plane();

  _index = -1;
  _channel = -1;
  _point1 = Point();
  _point2 = Point();

  time_slice = wire.GetTimeSlice();

  wire_all = wire.get_allwire();
}

MergeGeomWire::MergeGeomWire(int ident, const WireCell::MergeGeomWire& wire)
{
  _ident = ident;
  _plane = wire.plane();

  _index = -1;
  _channel = -1;
  _point1 = Point();
  _point2 = Point();

  time_slice = wire.GetTimeSlice();

  wire_all = wire.get_allwire();
}

MergeGeomWire::~MergeGeomWire(){
  wire_all.clear();
}

int MergeGeomWire::AddWire(const WireCell::GeomWire& wire){
  if (wire.plane() == _plane){
    auto it = find(wire_all.begin(),wire_all.end(),&wire);
    if (it == wire_all.end())
      wire_all.push_back(&wire);
    return 1;
  }else{
    return 0;
  }
}

int MergeGeomWire::AddWire(WireCell::MergeGeomWire& wire){
  if (wire.plane() == _plane){
 
    GeomWireSelection wires2 = wire.get_allwire();
    int flag = 0;
    const GeomWire *swire1, *swire2;
    
    for (int i=0;i!=wire_all.size();i++){
      swire1 = wire_all[i];

      for (int j=0;j!=wires2.size();j++){
	swire2 = wires2[j];
	
	//if(swire1->ident() == swire2->ident()){
	if(swire1->channel() == swire2->channel()){
	  flag = 1;
	  break;
	}
      }
      if (flag==1) break;
    }
    
    if (flag==1){
      for (int j=0;j!=wires2.size();j++){
	int flag1 = 0;
	swire2 = wires2[j];
	for (int i=0;i!=wire_all.size();i++){
	  swire1 = wire_all[i];
	  if (swire1->ident() == swire2->ident()){
	  //if (swire1->channel() == swire2->channel()){
	    flag1 = 1;
	    break;
	  }
	}
	if (flag1==0) wire_all.push_back(swire2);
      }
      
      return 1;
    }else{
      return 0;
    }

  }else{
    return 0;
  }
}
