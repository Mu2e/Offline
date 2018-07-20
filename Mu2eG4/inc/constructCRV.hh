#ifndef Mu2eG4_constructCRV_hh
#define Mu2eG4_constructCRV_hh
//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.hh,v 1.4 2012/11/19 23:03:24 genser Exp $
// $Author: genser $
// $Date: 2012/11/19 23:03:24 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructCRV(VolumeInfo   const & parent,
                    SimpleConfig const & _config
                    );

}

#endif /* Mu2eG4_constructCRV_hh */
