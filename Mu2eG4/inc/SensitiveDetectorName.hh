#ifndef SensitiveDetectorName_h
#define SensitiveDetectorName_h
// Define names of Sensitive Detectors
// 
// $Id: SensitiveDetectorName.hh,v 1.4 2011/03/09 19:49:32 genser Exp $
// $Author: genser $ 
// $Date: 2011/03/09 19:49:32 $
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

#endif
