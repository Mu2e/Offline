#ifndef constructDS_HH
#define constructDS_HH
//
// Free function to create the Detector Solenoid 
//
// $Id: constructDS.hh,v 1.1 2011/01/05 21:04:31 genser Exp $
// $Author: genser $
// $Date: 2011/01/05 21:04:31 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructDS(VolumeInfo   const & parent, 
                   SimpleConfig const * const _config
                   );

}

#endif
