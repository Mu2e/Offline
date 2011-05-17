#ifndef Mu2eG4_constructVirtualDetectors_hh
#define Mu2eG4_constructVirtualDetectors_hh
//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author KLG
//

namespace mu2e {

  class G4Helper;
  class SimpleConfig;

  void constructVirtualDetectors( SimpleConfig const * const _config
                                  );

}

#endif /* Mu2eG4_constructVirtualDetectors_hh */
