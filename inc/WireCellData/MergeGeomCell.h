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
    int GetTimeSlice(){ return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};

    bool GetContainTruthCell(){return contain_truth;};

    WireCell::GeomCellSelection get_allcell() const{ return cell_all;}
    WireCell::GeomCellSelection get_truthcell() const{return truth_cells;}

    bool CheckContainTruthCell(WireCell::CellChargeMap &ccmap);

  protected:
    WireCell::GeomCellSelection cell_all;
    
    int time_slice; // illustrate which time slice this is
    bool contain_truth; // whether it contain truth, default is not

    WireCell::GeomCellSelection truth_cells;
  };
}
#endif
