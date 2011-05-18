#ifndef Mu2eG4_SensitiveDetectorName_hh
#define Mu2eG4_SensitiveDetectorName_hh
// Define names of Sensitive Detectors
//
// $Id: SensitiveDetectorName.hh,v 1.6 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author KLG

namespace mu2e {

  class SensitiveDetectorName {

  public:

    // we define the functins here to avoid a run time undefined symbol error
    static char const * StrawGasVolume(){
      // string literals are statically allocated, so this is safe
      return "StrawGasVolume";
    }

    static char const * ItrackerGasVolume(){
       return "ItrackerGasVolume";
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

    static char const * StoppingTarget(){
      return "StoppingTarget";
    }

    static char const * CRSScintillatorBar(){
      return "CRSScintillatorBar";
    }

  };

} // namespace mu2e

#endif /* Mu2eG4_SensitiveDetectorName_hh */
