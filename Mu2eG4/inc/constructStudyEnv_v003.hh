#ifndef Mu2eG4_constructStudyEnv_v003_hh
#define Mu2eG4_constructStudyEnv_v003_hh
//
// Free function to create one of the study environment geometries
//
// $Id: constructStudyEnv_v003.hh,v 1.1 2013/04/05 19:52:40 genser Exp $
// $Author: genser $
// $Date: 2013/04/05 19:52:40 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructStudyEnv_v003(VolumeInfo   const & parentVInfo,
                              SimpleConfig const & _config
                              );

}

#endif /* Mu2eG4_constructStudyEnv_v003_hh */
