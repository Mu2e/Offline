#ifndef Mu2eG4_constructStudyEnv_v001_hh
#define Mu2eG4_constructStudyEnv_v001_hh
//
// Free function to create one of the study environment geometries
//
// $Id: constructStudyEnv_v001.hh,v 1.2 2013/04/05 19:51:34 genser Exp $
// $Author: genser $
// $Date: 2013/04/05 19:51:34 $
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
