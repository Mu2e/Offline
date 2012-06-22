#ifndef Mu2eG4_SensitiveDetectorName_hh
#define Mu2eG4_SensitiveDetectorName_hh
// Define names of Sensitive Detectors; revised to forward the names of the
// StepInstanceName names.
//
// $Id: SensitiveDetectorName.hh,v 1.13 2012/06/22 18:14:56 youzy Exp $
// $Author: youzy $
// $Date: 2012/06/22 18:14:56 $
//
// Original author KLG

#include "MCDataProducts/inc/StepInstanceName.hh"

namespace mu2e {

  class SensitiveDetectorName {

  public:

    // we define the functins here to avoid a run time undefined symbol error
    static char const * TrackerGas(){
      return StepInstanceName::name(StepInstanceName::tracker).c_str();
    }

    static char const * VirtualDetector(){
      return StepInstanceName::name(StepInstanceName::virtualdetector).c_str();
    }

    static char const * CaloCrystal(){
      return StepInstanceName::name(StepInstanceName::calorimeter).c_str();
    }

    static char const * CaloReadout(){
      return StepInstanceName::name(StepInstanceName::calorimeterRO).c_str();
    }

    static char const * ExtMonFNAL(){
      return StepInstanceName::name(StepInstanceName::ExtMonFNAL).c_str();
    }

    static char const * ExtMonUCITof(){
      return StepInstanceName::name(StepInstanceName::ExtMonUCITof).c_str();
    }

    static char const * StoppingTarget(){
      return StepInstanceName::name(StepInstanceName::stoppingtarget).c_str();
    }

    static char const * CRSScintillatorBar(){
      return StepInstanceName::name(StepInstanceName::CRV).c_str();
    }

    static char const * TTrackerDeviceSupport(){
      return StepInstanceName::name(StepInstanceName::ttrackerDS).c_str();
    }

    static char const * ProtonAbsorber() {
      return StepInstanceName::name(StepInstanceName::protonabsorber).c_str();
    }

    static char const * PSVacuum() {
      return StepInstanceName::name(StepInstanceName::PSVacuum).c_str();
    }

  };

} // namespace mu2e

#endif /* Mu2eG4_SensitiveDetectorName_hh */
