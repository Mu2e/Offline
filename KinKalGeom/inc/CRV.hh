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
      auto const& rightSector1() const { return *rightSector1_; }
      auto const& rightSector2() const { return *rightSector2_; }
      auto const& rightSector3() const { return *rightSector3_; }
      auto const& rightSector4() const { return *rightSector4_; }
      auto const& rightSector5() const { return *rightSector5_; }
      auto const& rightSector6() const { return *rightSector6_; }
      auto const& leftSector1() const { return *leftSector1_; }
      auto const& leftSector2() const { return *leftSector2_; }
      auto const& leftSector3() const { return *leftSector3_; }
      auto const& topSector1() const { return *topSector1_; }
      auto const& topSector2() const { return *topSector2_; }
      auto const& topSector3() const { return *topSector3_; }
      auto const& topSector4() const { return *topSector4_; }
      auto const& topSector5() const { return *topSector5_; }
      auto const& extSector1() const { return *extSector1_; }
      auto const& extSector2() const { return *extSector2_; }
      auto const& upstreamSector1() const { return *upstreamSector1_; }
      auto const& downstreamSector1() const { return *downstreamSector1_; }
      auto const& downstreamSector2() const { return *downstreamSector2_; }
      auto const& downstreamSector3() const { return *downstreamSector3_; }
      auto const& downstreamSector4() const { return *downstreamSector4_; }
      auto const& cryoSector1() const { return *cryoSector1_; }
      auto const& cryoSector2() const { return *cryoSector2_; }
      auto const& rightSector1Ptr() const { return rightSector1_; }
      auto const& rightSector2Ptr() const { return rightSector2_; }
      auto const& rightSector3Ptr() const { return rightSector3_; }
      auto const& rightSector4Ptr() const { return rightSector4_; }
      auto const& rightSector5Ptr() const { return rightSector5_; }
      auto const& rightSector6Ptr() const { return rightSector6_; }
      auto const& leftSector1Ptr() const { return leftSector1_; }
      auto const& leftSector2Ptr() const { return leftSector2_; }
      auto const& leftSector3Ptr() const { return leftSector3_; }
      auto const& topSector1Ptr() const { return topSector1_; }
      auto const& topSector2Ptr() const { return topSector2_; }
      auto const& topSector3Ptr() const { return topSector3_; }
      auto const& topSector4Ptr() const { return topSector4_; }
      auto const& topSector5Ptr() const { return topSector5_; }
      auto const& extSector1Ptr() const { return extSector1_; }
      auto const& extSector2Ptr() const { return extSector2_; }
      auto const& upstreamSector1Ptr() const { return upstreamSector1_; }
      auto const& downstreamSector1Ptr() const { return downstreamSector1_; }
      auto const& downstreamSector2Ptr() const { return downstreamSector2_; }
      auto const& downstreamSector3Ptr() const { return downstreamSector3_; }
      auto const& downstreamSector4Ptr() const { return downstreamSector4_; }
      auto const& cryoSector1Ptr() const { return cryoSector1_; }
      auto const& cryoSector2Ptr() const { return cryoSector2_; }
      private:
      RecPtr rightSector1_, rightSector2_, rightSector3_, rightSector4_, rightSector5_, rightSector6_, leftSector1_, leftSector2_, leftSector3_, topSector1_, topSector2_, topSector3_, topSector4_, topSector5_, extSector1_, extSector2_, upstreamSector1_, downstreamSector1_, downstreamSector2_, downstreamSector3_, downstreamSector4_, cryoSector1_, cryoSector2_;
    };
  }
}

#endif
