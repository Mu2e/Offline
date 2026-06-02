//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_CRV_hh
#define KinKalGeom_CRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include <vector>
#include <memory>
namespace mu2e {
  namespace KKGeom {
    using RecPtr = std::shared_ptr<KinKal::Rectangle>;
    struct KKCRVSector {
      std::string sname_;
      RecPtr sector_;
      double whw_;
    };
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
