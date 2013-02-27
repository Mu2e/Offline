#ifndef Mu2eG4_constructStudyEnv_v002_hh
#define Mu2eG4_constructStudyEnv_v002_hh
//
// Free function to create the Detector Solenoid
//
// $Id: constructStudyEnv_v002.hh,v 1.1 2013/02/27 03:48:50 genser Exp $
// $Author: genser $
// $Date: 2013/02/27 03:48:50 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructStudyEnv_v002(VolumeInfo   const & parentVInfo,
                              SimpleConfig const & _config
                              );

}

#endif /* Mu2eG4_constructStudyEnv_v002_hh */
