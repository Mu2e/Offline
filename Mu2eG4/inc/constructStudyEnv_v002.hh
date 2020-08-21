#ifndef Mu2eG4_constructStudyEnv_v002_hh
#define Mu2eG4_constructStudyEnv_v002_hh
//
// Free function to create one of the study environment geometries
//
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
