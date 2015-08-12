#ifndef WCVertex_h
#define WCVertex_h

#include "WireCellData/Point.h"
#include "WireCellData/MergeSpaceCell.h"
#include <vector>

namespace WireCell {
  class WCVertex {
    WCVertex(MergeSpaceCell& msc);
    ~WCVertex();
  public:
    Point Center();

  protected:
    Point center;
    MergeSpaceCell& msc;
  };
  
  typedef std::vector<WCVertex*> WCVertexSelection;
}

#endif
