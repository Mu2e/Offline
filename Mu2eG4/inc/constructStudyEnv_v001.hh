#ifndef Mu2eG4_constructStudyEnv_v001_hh
#define Mu2eG4_constructStudyEnv_v001_hh
//
// Free function to create the Detector Solenoid
//
// $Id: constructStudyEnv_v001.hh,v 1.1 2012/11/16 23:28:23 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:28:23 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructStudyEnv_v001(VolumeInfo   const & parentVInfo,
                              SimpleConfig const & _config
                              );

}

#endif /* Mu2eG4_constructStudyEnv_v001_hh */
