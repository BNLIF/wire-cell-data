#include "WCPData/TrackInfo.h"

using namespace WCP;

TrackInfo::TrackInfo(WCP::PointVector& tracking_path, std::vector<double>& dQ, std::vector<double>& dx, std::vector<double>& pu, std::vector<double>& pv, std::vector<double>& pw, std::vector<double>& pt, std::vector<double>& reduced_chi2)
  : tracking_path(tracking_path)
  , dQ(dQ)
  , dx(dx)
  , pu(pu)
  , pv(pv)
  , pw(pw)
  , pt(pt)
  , reduced_chi2(reduced_chi2)
  , parent_track(0)
{
  pcloud = new ToyPointCloud();
  pcloud->AddPoints(tracking_path);
  pcloud->build_kdtree_index();
}


TrackInfo::~TrackInfo(){
  delete pcloud;
}

std::pair<double, WCP::Point> TrackInfo::get_closest_point(WCP::Point &p){
  return pcloud->get_closest_point(p);
}

std::tuple<double, double, double> TrackInfo::get_closest_2d_dis(WCP::Point &p){
  std::pair<int, double> results_u = pcloud->get_closest_2d_dis(p, 0);
  std::pair<int, double> results_v = pcloud->get_closest_2d_dis(p, 1);
  std::pair<int, double> results_w = pcloud->get_closest_2d_dis(p, 2);

  return std::make_tuple(results_u.second, results_v.second, results_w.second);
}


bool TrackInfo::AddTrackInfo(WCP::TrackInfo* track, WCP::Point& pos, WCP::Point& self_pos, bool isparent){
  if (isparent){
    // parent tracks ...
    if (parent_track != 0) return false;
    parent_track = track;
    map_track_points[track] = std::make_pair(pos, self_pos);
  }else{
    //daughter tracks ...
    daughter_tracks.push_back(track);
    map_track_points[track] = std::make_pair(pos, self_pos);
  }
  return true;
}

std::pair<bool, std::pair<WCP::Point, WCP::Point> > TrackInfo::get_track_points(WCP::TrackInfo* track){
  auto it = map_track_points.find(track);

  if (it == map_track_points.end()){
    return std::make_pair(false,std::make_pair(WCP::Point(0,0,0), WCP::Point(0,0,0)) );
  }else{
    return std::make_pair(true, it->second);
  }
}
