#ifndef GeomWireCellData_WCPOINT_CLOUD_H
#define GeomWireCellData_WCPOINT_CLOUD_H

#include "WireCellData/SlimMergeGeomCell.h"

namespace WireCell{
  
  template <typename T>
    struct WCPointCloud
    {
      struct WCPoint
      {
	T  x,y,z;
	SlimMergeGeomCell *mcell;
	int index_u, index_v, index_w;
      };
      
      std::vector<WCPoint>  pts;
      
      // Must return the number of data points
      inline size_t kdtree_get_point_count() const { return pts.size(); }
      inline SlimMergeGeomCell* get_mcell(const size_t idx) const { return pts[idx].mcell;}
      
      // Returns the dim'th component of the idx'th point in the class:
      // Since this is inlined and the "dim" argument is typically an immediate value, the
      //  "if/else's" are actually solved at compile time.
      inline T kdtree_get_pt(const size_t idx, int dim) const
      {
	if (dim == 0) return pts[idx].x;
	else if (dim == 1) return pts[idx].y;
	else return pts[idx].z;
      }
      
      // Optional bounding-box computation: return false to default to a standard bbox computation loop.
      //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
      //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
      template <class BBOX>
      bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
      
    };
  
}

#endif
