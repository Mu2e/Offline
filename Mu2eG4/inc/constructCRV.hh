#ifndef Mu2eG4_constructCRV_hh
#define Mu2eG4_constructCRV_hh
//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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
