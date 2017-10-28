#ifndef WireCell_TOYPOINTCLOUD_H
#define WireCell_TOYPOINTCLOUD_H


#include "WireCellNanoflann/nanoflann.h" 
#include "WireCellData/WCPointCloud.h"
#include "WireCellQuickhull/QuickHull.h"
#include "WireCellQuickhull/MathUtils.h"

#include  <map>

typedef nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, WireCell::WCPointCloud<double> > ,
    WireCell::WCPointCloud<double>,
    3 /* dim */
  > my_kd_tree_t;

namespace WireCell{
  
  class ToyPointCloud {
  public:
    ToyPointCloud();
    ~ToyPointCloud();

    void AddPoint(WireCell::Point& p, std::tuple<int,int,int>& wires_index, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps, std::vector<std::tuple<int,int,int>>& wires_indices, WireCell::SlimMergeGeomCell *mcell);
    void build_kdtree_index();
    std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> get_hull();

    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, int N);
    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, double radius);
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, int N);
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, double radius);
    
    
    
    int get_num_points(){return cloud.pts.size();};
    std::vector<int>& get_mcell_indices(WireCell::SlimMergeGeomCell* mcell){
      return map_mcell_indices[mcell];
    };
    WireCell::WCPointCloud<double>& get_cloud(){return cloud;};

    void Print();
    
  protected:
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, double radius);

    
    WireCell::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
    // map of mcells to WCPoints
    std::map<WireCell::SlimMergeGeomCell*, std::vector<int>> map_mcell_indices; 
    
    

  };

  typedef std::vector<ToyPointCloud*> ToyPointCloudSelection;
   
}

#endif
