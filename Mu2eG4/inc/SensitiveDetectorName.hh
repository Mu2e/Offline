#ifndef Mu2eG4_SensitiveDetectorName_hh
#define Mu2eG4_SensitiveDetectorName_hh
// Define names of Sensitive Detectors
//
// $Id: SensitiveDetectorName.hh,v 1.11 2012/06/08 22:32:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/08 22:32:18 $
//
// Original author KLG

namespace mu2e {

  class SensitiveDetectorName {

  public:

    // we define the functins here to avoid a run time undefined symbol error
    static char const * TrackerGas(){
      // string literals are statically allocated, so this is safe
      return "tracker";
    }

    static char const * VirtualDetector(){
      return "VirtualDetector";
    }

    static char const * CaloCrystal(){
      return "CaloCrystal";
    }

    static char const * CaloReadout(){
      return "CaloReadout";
    }

    static char const * ExtMonFNAL(){
      return "ExtMonFNAL";
    }

    static char const * ExtMonUCITof(){
      return "ExtMonUCITof";
    }

    static char const * StoppingTarget(){
      return "StoppingTarget";
    }

    static char const * CRSScintillatorBar(){
      return "CRSScintillatorBar";
    }

    static char const * TTrackerDeviceSupport(){
      return "TTrackerDeviceSupport";
    }

    static char const * ProtonAbsorber() {
      return "ProtonAbsorber";
    }

  };

} // namespace mu2e

#endif /* Mu2eG4_SensitiveDetectorName_hh */
