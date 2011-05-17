#ifndef ConditionsService_unknownPDGIdName_hh
#define ConditionsService_unknownPDGIdName_hh
//
// Return a string identifying the name of a particle with an unknown PDGID.
//
// $Id: unknownPDGIdName.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
//

#include <string>

namespace mu2e {

  std::string unknownPDGIdName( unsigned int id);

}

#endif /* ConditionsService_unknownPDGIdName_hh */
