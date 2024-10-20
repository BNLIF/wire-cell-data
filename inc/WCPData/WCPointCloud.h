#ifndef GeomWCPData_WCPOINT_CLOUD_H
#define GeomWCPData_WCPOINT_CLOUD_H

#include "WCPData/SlimMergeGeomCell.h"

namespace WCP{
  
  template <typename T>
    struct WCPointCloud
    {
      struct WCPoint
      {
	T  x,y,z;
	SlimMergeGeomCell *mcell;
	int index_u, index_v, index_w;
	int index;
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


  template <typename T>
    struct WC2DPointCloud
    {
      struct WC2DPoint
      {
	T  x,y;
	int index;
      };
      
      std::vector<WC2DPoint>  pts;
      
      // Must return the number of data points
      inline size_t kdtree_get_point_count() const { return pts.size(); }
      
      // Returns the dim'th component of the idx'th point in the class:
      // Since this is inlined and the "dim" argument is typically an immediate value, the
      //  "if/else's" are actually solved at compile time.
      inline T kdtree_get_pt(const size_t idx, int dim) const
      {
	      if (dim == 0) return pts[idx].x;
	      else if (dim == 1) return pts[idx].y;
        return pts[idx].x;
      }
      
      // Optional bounding-box computation: return false to default to a standard bbox computation loop.
      //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
      //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
      template <class BBOX>
      bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
      
    };


  template <typename T>
    struct CTPointCloud
    {
      struct CTPoint
      {
	T  x,y;
	int channel;
	int time_slice;
	int charge;
	int charge_err;
	int index;
      };
      
      std::vector<CTPoint>  pts;
      
      // Must return the number of data points
      inline size_t kdtree_get_point_count() const { return pts.size(); }
      
      // Returns the dim'th component of the idx'th point in the class:
      // Since this is inlined and the "dim" argument is typically an immediate value, the
      //  "if/else's" are actually solved at compile time.
      inline T kdtree_get_pt(const size_t idx, int dim) const
      {
      	if (dim == 0) return pts[idx].x;
	      else if (dim == 1) return pts[idx].y;
        return pts[idx].x;
      }
      
      // Optional bounding-box computation: return false to default to a standard bbox computation loop.
      //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
      //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
      template <class BBOX>
      bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
      
    };

  
}

#endif
