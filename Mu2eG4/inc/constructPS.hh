#ifndef Mu2eG4_constructPS_hh
#define Mu2eG4_constructPS_hh
//
// Free function to create  Production Solenoid and Production Target
//
// $Id: constructPS.hh,v 1.5 2012/03/05 19:38:17 genser Exp $
// $Author: genser $
// $Date: 2012/03/05 19:38:17 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;
  class SensitiveDetectorHelper;

  void constructPS( const VolumeInfo& parent,
                    const SimpleConfig& _config,
                    const SensitiveDetectorHelper& sdHelper
                    );

}

#endif /* Mu2eG4_constructPS_hh */
