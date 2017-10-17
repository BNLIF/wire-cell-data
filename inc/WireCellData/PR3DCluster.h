#include "WireCellData/SlimMergeGeomCell.h"



namespace WireCell{
  class PR3DCluster{
  public:
    PR3DCluster(int cluster_id);
    ~PR3DCluster();
    void AddCell(SlimMergeGeomCell* mcell, int time_slice);
    // void AddCell(SlimMergeGeomCell* mcell, int *time_slices, int ntime_slice);
    int get_cluster_id(){return cluster_id;};
    int get_num_mcells(){return mcells.size();};
    int get_num_time_slices(){return time_cells_set_map.size();};

    void Remove_duplicated_mcells();
    SMGCSelection Is_Connected(PR3DCluster* cluster1, int offset);
    std::map<int,SMGCSet>& get_time_cells_set_map(){return time_cells_set_map;};
    
    
  protected:
    int cluster_id;
    SMGCSelection mcells;
    /* std::map<int,SMGCSelection> time_cells_map; */
    /* std::map<SlimMergeGeomCell*, std::vector<int>> cell_times_map; */
    
    std::map<int,SMGCSet> time_cells_set_map;
    std::map<SlimMergeGeomCell*, std::set<int>> cell_times_set_map;   
  };
  typedef std::vector<PR3DCluster*> PR3DClusterSelection;
}
