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

    void AddPoint(WireCell::Point& p, WireCell::SlimMergeGeomCell *mcell);
    void AddPoints(WireCell::PointVector& ps, WireCell::SlimMergeGeomCell *mcell);
    void build_kdtree_index();
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, int N);
    std::vector<std::pair<size_t,double>> get_closest_index(WireCell::Point& p, double radius);

    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, int N);
    std::map<WireCell::SlimMergeGeomCell*, WireCell::Point> get_closest_mcell(WireCell::Point& p, double radius);
    
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, int N);
    std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>> get_closest_points(WireCell::Point& p, double radius);
    
    std::vector<std::pair<WireCell::SlimMergeGeomCell*, WireCell::Point>> get_hull();
    int get_num_points(){return cloud.pts.size();};

    int get_wcpoint_index(WCPointCloud<double>::WCPoint* wcpoint){return map_wcpoint_index[wcpoint];};
    std::vector<WCPointCloud<double>::WCPoint*>& get_mcell_wcpoints(WireCell::SlimMergeGeomCell* mcell){return map_mcell_wcpoints[mcell];};
    
  protected:
    WireCell::WCPointCloud<double> cloud;
    my_kd_tree_t *index;
    
    // map of mcells to WCPoints
    std::map<WireCell::SlimMergeGeomCell*, std::vector<WCPointCloud<double>::WCPoint*>> map_mcell_wcpoints; 
    // map of WCPoint to index in the cloud 
    std::map<WCPointCloud<double>::WCPoint*,int> map_wcpoint_index;
    

  };

  typedef std::vector<ToyPointCloud*> ToyPointCloudSelection;
   
}

#endif
