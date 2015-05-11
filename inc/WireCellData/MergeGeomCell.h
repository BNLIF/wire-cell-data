#ifndef GeomWireCellData_MergeCell_h
#define GeomWireCellData_MergeCell_h

#include "WireCellData/GeomCell.h"
namespace WireCell{
  class MergeGeomCell : public WireCell::GeomCell {
  public: 
    MergeGeomCell(int ident, const WireCell::GeomCell& cell);
    MergeGeomCell(int ident, const WireCell::MergeGeomCell& cell);
    ~MergeGeomCell();

    double cross_section() const;
    Point center() const;

    int AddCell(const WireCell::GeomCell& cell);
    int AddCell(WireCell::MergeGeomCell& cell);

    WireCell::GeomCellSelection get_allcell() const{ return cell_all;}

    protected:
    WireCell::GeomCellSelection cell_all;
  };
}
#endif
