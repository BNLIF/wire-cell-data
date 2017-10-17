#include "WireCellData/SlimMergeGeomCell.h"



namespace WireCell{
  class PR3DCluster{
  public:
    PR3DCluster(int cluster_id);
    ~PR3DCluster();
  protected:
    int cluster_id;
    SMGCSelection mcells;
    std::map<int,SMGCSelection> time_cells_map;
    std::map<SlimMergeGeomCell*, std::vector<int>> cell_times_map;
  };
}
