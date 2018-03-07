#ifndef Mu2eG4_constructStoppingTarget_hh
#define Mu2eG4_constructStoppingTarget_hh
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.5 2011/08/04 18:52:25 genser Exp $
// $Author: genser $
// $Date: 2011/08/04 18:52:25 $
//
// Original author Peter Shanahan
//
// Notes:

namespace mu2e{

  class VolumeInfo;
  class SimpleConfig;
  class SensitiveDetectorHelper;

  VolumeInfo constructStoppingTarget(const VolumeInfo& parent,
                                     const SimpleConfig& _config,
                                     const SensitiveDetectorHelper& sdHelper
                                     );

}  // end namespace mu2e

#endif /* Mu2eG4_constructStoppingTarget_hh */
