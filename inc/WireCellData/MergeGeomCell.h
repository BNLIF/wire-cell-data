#ifndef GeomWireCellData_MergeCell_h
#define GeomWireCellData_MergeCell_h

#include "WireCellData/GeomCell.h"
namespace WireCell{
  class MergeGeomCell : public WireCell::GeomCell {
  public: 
    MergeGeomCell(int ident, const WireCell::GeomCell cell);
    
    int AddCell(const WireCell::GeomCell cell);

    WireCell::GeomCellSelection get_allcell(){ return cell_all;}
  protected:
    WireCell::GeomCellSelection cell_all;
  };
}
#endif
