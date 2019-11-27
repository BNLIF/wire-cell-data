#ifndef WCPTrackInfo_h
#define WCPTrackInfo_h

#include "WCPData/Point.h"
#include "WCPData/ToyPointCloud.h"
#include <vector>
#include <map>

namespace WCP{

  class TrackInfo {
  public:
    TrackInfo(WCP::PointVector& tracking_path, std::vector<double>& dQ, std::vector<double>& dx, std::vector<double>& pu, std::vector<double>& pv, std::vector<double>& pw, std::vector<double>& pt, std::vector<double>& reduced_chi2);
    ~TrackInfo();

    WCP::PointVector& get_tracking_path(){return tracking_path;};
    std::vector<double>& get_dQ(){return dQ;};
    std::vector<double>& get_dx(){return dx;};
    std::vector<double>& get_pu(){return pu;};
    std::vector<double>& get_pv(){return pv;};
    std::vector<double>& get_pw(){return pw;};
    std::vector<double>& get_pt(){return pt;};
    std::vector<double>& get_reduced_chi2(){return reduced_chi2;};

    std::pair<double, WCP::Point> get_closest_point(WCP::Point &p); 
    
    WCP::TrackInfo* get_parent_track(){return parent_track;};
    std::vector<WCP::TrackInfo*> get_daughter_tracks(){return daughter_tracks;};

   bool AddTrackInfo(WCP::TrackInfo* track, WCP::Point& pos, WCP::Point& self_pos, bool isparent=false);

   std::pair<bool, std::pair<WCP::Point, WCP::Point> > get_track_points(WCP::TrackInfo* track);
   
  protected:
    WCP::PointVector tracking_path;
    std::vector<double> dQ;
    std::vector<double> dx;
    std::vector<double> pu;
    std::vector<double> pv;
    std::vector<double> pw;
    std::vector<double> pt;
    std::vector<double> reduced_chi2;

    // point cloud ...
    WCP::ToyPointCloud* pcloud;
    
    // parent information
    WCP::TrackInfo* parent_track;
    std::vector<WCP::TrackInfo*> daughter_tracks;

    // 
    std::map<WCP::TrackInfo*, std::pair<WCP::Point, WCP::Point> > map_track_points;
    
  };

  typedef std::vector<WCP::TrackInfo*> TrackInfoSelection;
  
}

#endif
