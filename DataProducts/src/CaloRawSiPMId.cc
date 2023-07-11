#include "Offline/DataProducts/inc/CaloRawSiPMId.hh"
#include <iostream>

namespace mu2e {

  std::ostream& operator<<(std::ostream& ost,
                           const CaloRawSiPMId& id ){
    ost << id.id() ;
    return ost;
  }


} // end namespace mu2e
