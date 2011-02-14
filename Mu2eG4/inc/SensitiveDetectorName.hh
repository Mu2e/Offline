#ifndef SensitiveDetectorName_h
#define SensitiveDetectorName_h
// Define names of Sensitive Detectors
// 
// $Id: SensitiveDetectorName.hh,v 1.2 2011/02/14 23:20:01 logash Exp $
// $Author: logash $ 
// $Date: 2011/02/14 23:20:01 $
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

    static char const * StoppingTarget(){
      return "StoppingTarget";
    }

  };

} // namespace mu2e

#endif
