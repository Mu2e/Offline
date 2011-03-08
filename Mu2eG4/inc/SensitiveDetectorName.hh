#ifndef SensitiveDetectorName_h
#define SensitiveDetectorName_h
// Define names of Sensitive Detectors
// 
// $Id: SensitiveDetectorName.hh,v 1.3 2011/03/08 08:31:51 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2011/03/08 08:31:51 $
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

  };

} // namespace mu2e

#endif
