#ifndef Mu2eG4_constructProtonAbsorber_hh
#define Mu2eG4_constructProtonAbsorber_hh
//
// Free function to construct Proton Absorber
//
// $Id: constructProtonAbsorber.hh,v 1.3 2012/11/19 23:03:24 genser Exp $
// $Author: genser $
// $Date: 2012/11/19 23:03:24 $
//
// Original author KLG
//

namespace mu2e {

  class SimpleConfig;
  class SensitiveDetectorHelper;

  void constructProtonAbsorber(const SimpleConfig& _config,
                               const SensitiveDetectorHelper& sdHelper
                               );

}

#endif /* Mu2eG4_constructProtonAbsorber_hh */
