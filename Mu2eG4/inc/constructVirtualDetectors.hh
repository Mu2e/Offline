#ifndef Mu2eG4_constructVirtualDetectors_hh
#define Mu2eG4_constructVirtualDetectors_hh
//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.hh,v 1.3 2012/06/05 16:19:24 genser Exp $
// $Author: genser $
// $Date: 2012/06/05 16:19:24 $
//
// Original author KLG
//

namespace mu2e {

  class SimpleConfig;
  class SensitiveDetectorHelper;

  void constructVirtualDetectors( const SimpleConfig& _config,
                                  const SensitiveDetectorHelper& sdHelper
                                  );

}

#endif /* Mu2eG4_constructVirtualDetectors_hh */
