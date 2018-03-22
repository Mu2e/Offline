#ifndef Mu2eG4_constructDS_hh
#define Mu2eG4_constructDS_hh
//
// Free function to create the Detector Solenoid
//
// $Id: constructDS.hh,v 1.4 2012/11/16 23:45:24 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:45:24 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructDS(const VolumeInfo& parent,
                   const SimpleConfig& config
                   );

}

#endif /* Mu2eG4_constructDS_hh */
