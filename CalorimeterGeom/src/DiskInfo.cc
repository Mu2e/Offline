#include "Offline/CalorimeterGeom/inc/DiskInfo.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  DiskInfo::DiskInfo() { recompute_(); }

  CLHEP::Hep3Vector DiskInfo::toGlobal(const CLHEP::Hep3Vector& local)  const
  {
    return origin_ + inverseRotation_*local;
  }

  CLHEP::Hep3Vector DiskInfo::toLocal (const CLHEP::Hep3Vector& global) const
  {
    return rotation_*(global - origin_);
  }

  CLHEP::Hep3Vector DiskInfo::toLocalFF (const CLHEP::Hep3Vector& g) const
  {
    return toLocal(g) - originToCrystalOrigin_;
  }

  CLHEP::Hep3Vector DiskInfo::toGlobalFF(const CLHEP::Hep3Vector& l) const
  {
    return toGlobal(l + originToCrystalOrigin_);
  }

  void DiskInfo::recompute_() {
    inverseRotation_  = rotation_.inverse();
    frontFaceCenter_  = origin_ + inverseRotation_*ffLocal_;
    backFaceCenter_   = origin_ + inverseRotation_*bfLocal_;
    crystalDirection_ = inverseRotation_*dirLocal_;
  }

  void DiskInfo::print(std::ostream &os) const {
     os<<"origin Mu2e           "<<origin_<<std::endl;
     os<<"size                  "<<size_<<std::endl;
     os<<"rotation              "<<rotation_<<std::endl;
     os<<"Front face            "<<frontFaceCenter_<<std::endl;
     os<<"Back face             "<<backFaceCenter_<<std::endl;
     os<<"Crystal direction     "<<crystalDirection_<<std::endl;
     os<<"origin local          "<<originLocal_<<std::endl;
     os<<"originToCrystalOrigin "<<originToCrystalOrigin_<<std::endl;
     os<<"Front face local      "<<ffLocal_<<std::endl;
     os<<"Back face local       "<<bfLocal_<<std::endl;
     os<<"direction local       "<<dirLocal_<<std::endl;
     os<<"FEBZOffset            "<<FEBZOffset_<<std::endl;
     os<<"FEBZLength            "<<FEBZLength_<<std::endl;
  }


}
