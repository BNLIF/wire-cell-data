#ifndef GeomWireCellData_MergeCell_h
#define GeomWireCellData_MergeCell_h

#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWireCellMap.h"

namespace WireCell{
  class MergeGeomCell : public WireCell::GeomCell {
  public: 
    MergeGeomCell();
    MergeGeomCell(int ident, const WireCell::GeomCell& cell);
    MergeGeomCell(int ident, const WireCell::MergeGeomCell& cell);
    ~MergeGeomCell();

    double cross_section() const;
    Point center() const;

    bool Connected(const WireCell::GeomCell& cell1,const WireCell::GeomCell& cell2);
    void AddNewCell(const WireCell::GeomCell& cell);
    
    int AddCell(const WireCell::GeomCell& cell, double dis = 0.0002 * units::m);
    int AddCell(WireCell::MergeGeomCell& cell, double dis = 0.0002 * units::m);


    int GetTimeSlice() const { return time_slice;};
    void SetTimeSlice(int time){time_slice = time;};

    bool Overlap(const WireCell::MergeGeomCell& cell, float num = 0.1) const;
    int  Overlap1(const WireCell::MergeGeomCell& cell, float num = 0.1) const;
    bool GetContainTruthCell()const {return contain_truth;};
    float GetTruthCharge() const {return truth_charge;};
    void FindEdges();

    //function to find out corner cells, 
    void FindCorners(WireCell::GeomCellMap& cmap, WireCell::GeomWireMap& wmap);

    WireCell::GeomCellSelection get_allcell() const { return cell_all;}
    WireCell::GeomCellSelection get_edgecells() const { return edge_cells;}
    WireCell::GeomCellSelection get_truthcell() const{return truth_cells;}
    WireCell::GeomCellSelection get_cornercells()const { return corner_cells;}
    //std::vector<int>& get_cornercells_index() {return corner_cells_index;}
    WireCell::GeomCellSelection get_cornercells(int index1,int index2)const { return corner_cells_group[index1][index2];}
    CellIndexMap get_cornercells_index(){return corner_cells_index;}
      
    void Organize_edge_boundary();

    bool CheckContainTruthCell(WireCell::CellChargeMap &ccmap);
    bool IsBlob(){return blob;};
    bool IsSimpleBlob(){return simple_blob;};
    void SetSimpleBlob(bool value){simple_blob = value;};

    int index(int index1, int index2);
    int index1(int index);
    int index2(int index);

    void set_id(int num) const {mcell_id = num;};
    int get_id() const{return mcell_id;};

    GeomWireSelection get_uwires() const{return uwires;};
    GeomWireSelection get_vwires() const{return vwires;};
    GeomWireSelection get_wwires() const{return wwires;};
    GeomWireSelection get_edge_wires() const{return edge_wires;};
    GeomCellSelection& get_edge_cells(){return edge_cells;};

    //hack for now
    std::set<int> ewires;
    
  protected:
    mutable int mcell_id;

    WireCell::GeomWireSelection uwires;
    WireCell::GeomWireSelection vwires;
    WireCell::GeomWireSelection wwires;


    WireCell::GeomCellSelection cell_all;
    WireCell::GeomCellSelection edge_cells;
    WireCell::GeomWireSelection edge_wires;

    WireCell::GeomCellSelection corner_cells;
    WireCell::GeomCellSelection corner_cells_group[6][6];

    //std::vector<int> corner_cells_index;
    CellIndexMap corner_cells_index;
    
    int time_slice; // illustrate which time slice this is
    bool contain_truth; // whether it contain truth, default is not
    float truth_charge;
    bool blob;
    bool simple_blob;
    bool flag_corner;
    bool flag_edge;

    mutable int flag_center;
    mutable int flag_cross_section;
    mutable double area;
    mutable Point ret;

    WireCell::GeomCellSelection truth_cells;
  };

   /// Compare ident
    struct MergeGeomCellCompare {
      bool operator() (const MergeGeomCell* a, const MergeGeomCell* b) const {
	if (a->cross_section() == b->cross_section()) {
	  return a<b;
	}
	return a->cross_section() < b->cross_section();
      }
      
    };

    /// Used to store a definitive, ordered set of cells
    typedef std::set<WireCell::MergeGeomCell*, MergeGeomCellCompare> MergeGeomCellSet;
}
#endif
