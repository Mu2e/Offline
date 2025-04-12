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
  namespace KinKalGeom {
    class CRV {
    public:
      using RecPtr = std::shared_ptr<KinKal::Rectangle>;
      // default constructor with nominal geometry
      CRV();
      // accessors
      // return by reference
      auto const& rightSectorBlock() const { return *rightSectorBlock_; }
      auto const& rightSectorHole() const { return *rightSectorHole_; }
      auto const& cryoSectorHole() const { return *cryoSectorHole_; }
      auto const& leftSectorBlock() const { return *leftSectorBlock_; }
      auto const& leftSectorHole() const { return *leftSectorHole_; }
      auto const& topSectorBlock() const { return *topSectorBlock_; }
      auto const& extSectorBlock() const { return *extSectorBlock_; }
      auto const& upstreamSector1() const { return *upstreamSectorBlock_; }
      auto const& downstreamSectorBlock() const { return *downstreamSectorBlock_; }
      auto const& downstreamSectorHole() const { return *downstreamSectorHole_; }
      auto const& cryoSector1() const { return *cryoSector1_; }
      auto const& cryoSector2() const { return *cryoSector2_; }
      //
      auto const& rightSectorBlockPtr() const { return rightSectorBlock_; }
      auto const& rightSectorHolePtr() const { return rightSectorHole_; }
      auto const& cryoSectorHolePtr() const { return cryoSectorHole_; }
      auto const& leftSectorBlockPtr() const { return leftSectorBlock_; }
      auto const& leftSectorHolePtr() const { return leftSectorHole_; }
      auto const& topSectorBlockPtr() const { return topSectorBlock_; }
      auto const& extSectorBlockPtr() const { return extSectorBlock_; }
      auto const& upstreamSectorBlockPtr() const { return upstreamSectorBlock_; }
      auto const& downstreamSectorBlockPtr() const { return downstreamSectorBlock_; }
      auto const& downstreamSectorHolePtr() const { return downstreamSectorHole_; }
      auto const& cryoSector1Ptr() const { return cryoSector1_; }
      auto const& cryoSector2Ptr() const { return cryoSector2_; }
      private:
      RecPtr rightSectorBlock_, rightSectorHole_, cryoSectorHole_, leftSectorBlock_, leftSectorHole_, topSectorBlock_, extSectorBlock_, upstreamSectorBlock_, downstreamSectorBlock_, downstreamSectorHole_, cryoSector1_, cryoSector2_;
    };
  }
}

#endif
