#ifndef WCVertex_h
#define WCVertex_h
#include "WireCellData/Point.h"

namespace WireCell {
  class WCVertex {
  public:
    Point Center();

  protected:
    Point center;
  };
}

#endif
