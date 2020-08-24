#ifndef WIRECELL_PHOTON_LIBRARY_H
#define WIRECELL_PHOTON_LIBRARY_H

#include <vector>
#include <map>
#include <list>

namespace WCP{
  class Photon_Library {
  public:
    std::map<int,int> map_lib_pmt, map_pmt_lib;
    std::vector<std::list<std::pair<int,float>>> library;
    double scaling_light_mag, rel_light_yield_err;
    
    Photon_Library(double eventTime, int run_no = 0, bool flag_data = true, bool flag_add_light_yield_err = false, bool flag_timestamp = false);
  };
  
}

#endif
