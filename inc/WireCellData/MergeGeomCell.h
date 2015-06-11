#ifndef GeomWireCellData_MergeCell_h
#define GeomWireCellData_MergeCell_h

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWireCellMap.h"

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
    int GetTimeSlice() const { return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};

    bool Overlap(const WireCell::MergeGeomCell& cell) const;
    bool GetContainTruthCell()const {return contain_truth;};
    float GetTruthCharge() const {return truth_charge;};
    void FindEdges();

    //function to find out corner cells, 
    void FindCorners(WireCell::GeomCellMap& cmap, WireCell::GeomWireMap& wmap);

    WireCell::GeomCellSelection get_allcell() const{ return cell_all;}
    WireCell::GeomCellSelection get_edgecells() const{ return edge_cells;}
    WireCell::GeomCellSelection get_truthcell() const{return truth_cells;}
    WireCell::GeomCellSelection get_cornercells()const { return corner_cells;}
    std::vector<int>& get_cornercells_index() {return corner_cells_index;}

    bool CheckContainTruthCell(WireCell::CellChargeMap &ccmap);

  protected:
    WireCell::GeomCellSelection cell_all;
    WireCell::GeomCellSelection edge_cells;

    WireCell::GeomWireSelection edge_wires;
    WireCell::GeomCellSelection corner_cells;
    std::vector<int> corner_cells_index;
    
    int time_slice; // illustrate which time slice this is
    bool contain_truth; // whether it contain truth, default is not
    float truth_charge;

    WireCell::GeomCellSelection truth_cells;
  };

   /// Compare ident
    struct MergeGeomCellCompare {
      bool operator() (const MergeGeomCell* a, const MergeGeomCell* b) const {
	return a->cross_section() < b->cross_section();
      }
      
    };

    /// Used to store a definitive, ordered set of cells
    typedef std::set<WireCell::MergeGeomCell*, MergeGeomCellCompare> MergeGeomCellSet;
}
#endif
