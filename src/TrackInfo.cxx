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
