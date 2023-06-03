#include "Offline/DataProducts/inc/CaloId.hh"

namespace mu2e {


  uint16_t CaloId::disk() {
      if(crystal()<_nCrystalPerDisk ||
         (_id >= _nCrystalChannel && _id < _nCrystalChannel+_nPINDiodPerDisk) ) {
        return 0;
      } else {
        return 1;
      }
    }

  bool CaloId::isCaphri() {
      if(std::find(_caphriId.begin(),_caphriId.end(),crystal())==_caphriId.end()) return false;
      return true;
    }

} // end namespace mu2e
