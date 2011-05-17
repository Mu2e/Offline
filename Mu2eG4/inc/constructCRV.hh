#ifndef Mu2eG4_constructCRV_hh
#define Mu2eG4_constructCRV_hh
//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructCRV(VolumeInfo   const & parent, 
                    SimpleConfig const * const _config
                    );

}

#endif /* Mu2eG4_constructCRV_hh */
