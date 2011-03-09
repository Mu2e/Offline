#ifndef constructCRV_HH
#define constructCRV_HH
//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
// $Id: constructCRV.hh,v 1.1 2011/03/09 19:25:47 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:25:47 $
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

#endif
