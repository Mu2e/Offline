#ifndef constructVirtualDetectors_HH
#define constructVirtualDetectors_HH
//
// Free function to create the virtual detectors
//
// $Id: constructVirtualDetectors.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

namespace mu2e {

  class G4Helper;
  class SimpleConfig;

  void constructVirtualDetectors( SimpleConfig const * const _config
                                  );

}

#endif
