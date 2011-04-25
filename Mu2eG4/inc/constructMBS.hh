#ifndef constructMBS_HH
#define constructMBS_HH
//
// Free function to create Muon Beam Stop and some elements of the Cryostat in G4
//
// $Id: constructMBS.hh,v 1.1 2011/04/25 19:17:07 genser Exp $
// $Author: genser $
// $Date: 2011/04/25 19:17:07 $
//
// Original author KLG
//

namespace mu2e {
  
  class SimpleConfig;

  void constructMBS(SimpleConfig const * const _config);

}

#endif
