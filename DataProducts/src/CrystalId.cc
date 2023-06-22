#include "Offline/DataProducts/inc/CrystalId.hh"
#include <algorithm>
#include <iostream>

namespace mu2e {

  bool CrystalId::isCaphri() const {
    if(std::find(CaloConst::_caphriId.begin(),
                 CaloConst::_caphriId.end(),_id)
       == CaloConst::_caphriId.end()) {
      return false;
    } else {
      return true;
    }
  }

  std::ostream& operator<<(std::ostream& ost,
                           const CrystalId& id ){
    ost << id.id() ;
    return ost;
  }


} // end namespace mu2e
