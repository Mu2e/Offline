#ifndef Mu2eG4_customizeChargedPionDecay_hh
#define Mu2eG4_customizeChargedPionDecay_hh
//
// Customize the e nu decay channel of charged pions.
//
// PDG - Set the branching fraction to the PDG value
// Off - Set the branching fraction to 0.
// All - Set the branching fraction to 100%
// A numerical value - Set the branching fraction to the given number.
//
// $Id: customizeChargedPionDecay.hh,v 1.1 2012/07/10 21:16:53 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/10 21:16:53 $
//

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e{

  void customizeChargedPionDecay(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

}  // end namespace mu2e

#endif /* Mu2eG4_customizeChargedPionDecay_hh */
