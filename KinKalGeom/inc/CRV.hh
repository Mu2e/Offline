//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_CRV_hh
#define KinKalGeom_CRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/General/Vectors.hh"
#include <cmath>
#include <string>
#include <utility>
#include <vector>
#include <memory>
namespace mu2e {
  namespace KKGeom {
    using RecPtr = std::shared_ptr<KinKal::Rectangle>;
    // one scintillator layer of a sector: its offset from the sector midplane and the spans its bars cover across the bars
    struct KKCRVLayer {
      double w_ = 0.0; // layer center along the sector normal, relative to the sector midplane (mm)
      std::vector<std::pair<double,double>> v_; // [low,high] along the sector vDirection, relative to the sector center (mm), one per module
    };
    struct KKCRVSector {
      std::string sname_;
      RecPtr sector_; // midplane, spanning every layer: adjacent sectors overlap where their staggered layers interleave
      double whw_;
      std::vector<KKCRVLayer> layers_;
      // fraction of the layers that a straight line through pos (on the midplane) along dir crosses; 1 if no layers are given
      double layerFraction(KinKal::VEC3 const& pos, KinKal::VEC3 const& dir) const;
    };
    inline double KKCRVSector::layerFraction(KinKal::VEC3 const& pos, KinKal::VEC3 const& dir) const {
      if(layers_.empty()) return 1.0;
      double const dw = dir.Dot(sector_->normal());
      double const slope = std::fabs(dw) > 1.0e-6 ? dir.Dot(sector_->vDirection())/dw : 0.0; // dv/dw along the line
      double const v0 = (pos - sector_->center()).Dot(sector_->vDirection());
      unsigned ncross = 0;
      for(auto const& layer : layers_){
        double const v = v0 + slope*layer.w_;
        for(auto const& span : layer.v_){
          if(v >= span.first && v <= span.second){ ++ncross; break; }
        }
      }
      return ncross/double(layers_.size());
    }
    class CRV {
      public:
        // construct with rectangles representing the midplane of CRV sectors (shields)
        CRV(std::vector<KKCRVSector> const& sectors) : sectors_(sectors) {}
        // accessors
        auto const& sectors() const { return sectors_; }
        auto const& sector(size_t isect) const { return sectors_[isect].sector_; }
        auto const& sectorName(size_t isect) const { return sectors_[isect].sname_; }
        auto const& sectorHalfWidth(size_t isect) const { return sectors_[isect].whw_; }
      private:
        std::vector<KKCRVSector> sectors_;
    };
  }
}

#endif
