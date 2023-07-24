#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include <iostream>

namespace mu2e {

  CaloSiPMId::value_type CaloSiPMId::disk() const {
    if(_id<CaloConst::_nCrystalPerDisk*CaloConst::_nSiPMPerCrystal ||
         (_id >= CaloConst::_nCrystalChannel &&
          _id < CaloConst::_nCrystalChannel+CaloConst::_nPINDiodPerDisk) ) {
        return 0;
      } else {
        return 1;
      }
    }

  CaloSiPMId::value_type CaloSiPMId::detType() const {
    if(crystal().isCaphri()) {
      return CaloConst::detType::CAPHRI;
    } else if(isCrystal()) {
      return CaloConst::detType::CsI;
    } else {
      return CaloConst::detType::PINDiode;
    }
  }

  std::ostream& operator<<(std::ostream& ost,
                           const CaloSiPMId& id ){
    ost << id.id() ;
    return ost;
  }


} // end namespace mu2e
