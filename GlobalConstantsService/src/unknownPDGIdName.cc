//
// Return a string identifying the name of a particle with an unknown PDGID.
//
// Original author Rob Kutschke
//

#include <sstream>

#include "Offline/GlobalConstantsService/inc/unknownPDGIdName.hh"

namespace mu2e {

  std::string unknownPDGIdName( unsigned int id ){
    std::ostringstream os;
    os << "UnknownPDGId_" << id;
    return os.str();
  }

}
