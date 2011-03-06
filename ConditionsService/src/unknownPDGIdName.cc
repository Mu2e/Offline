//
// Return a string identifying the name of a particle with an unknown PDGID.
//
// $Id: unknownPDGIdName.cc,v 1.1 2011/03/06 00:31:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/03/06 00:31:43 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "ConditionsService/inc/unknownPDGIdName.hh"

namespace mu2e {

  std::string unknownPDGIdName( unsigned int id ){
    std::ostringstream os;
    os << "UnknownPDGId_" << id;
    return os.str();
  }

}
