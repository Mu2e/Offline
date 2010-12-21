#ifndef SensitiveDetectorName_h
#define SensitiveDetectorName_h
// Define names of Sensitive Detectors
// 
// $Id: SensitiveDetectorName.hh,v 1.1 2010/12/21 21:33:47 genser Exp $
// $Author: genser $ 
// $Date: 2010/12/21 21:33:47 $
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

    static char const * VirtualDetector(){
      return "VirtualDetector";
    }

    static char const * CaloCrystal(){
      return "CaloCrystal";
    }

    static char const * CaloReadout(){
      return "CaloReadout";
    }

  };

} // namespace mu2e

#endif
