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

  CaloConst::detType CaloSiPMId::detType() const {
    if(isCrystal()){
      if(crystal().isCaphri()) {
        return CaloConst::detType::CAPHRI;
      } else {
        return CaloConst::detType::CsI;
      }
    } else if(isPINDiode()) {
      return CaloConst::detType::PINDiode;
    } else {
      return CaloConst::detType::Invalid;
    }
  }

  std::ostream& operator<<(std::ostream& ost,
                           const CaloSiPMId& id ){
    ost << id.id() ;
    return ost;
  }


} // end namespace mu2e
