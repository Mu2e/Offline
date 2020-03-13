#ifndef Mu2eG4_SensitiveDetectorName_hh
#define Mu2eG4_SensitiveDetectorName_hh
// Define names of Sensitive Detectors; revised to forward the names of the
// StepInstanceName names.
//
// $Id: SensitiveDetectorName.hh,v 1.17 2014/06/05 21:06:34 genser Exp $
// $Author: genser $
// $Date: 2014/06/05 21:06:34 $
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

    static char const * CaloReadoutCard(){
      return StepInstanceName::name(StepInstanceName::calorimeterROCard).c_str();
    }

    static char const * CaloCrate(){
      return StepInstanceName::name(StepInstanceName::calorimeterCrate).c_str();
    }

    static char const * ExtMonFNAL(){
      return "ExtMonFNAL";
    }

    static char const * ExtMonUCITof(){
      return StepInstanceName::name(StepInstanceName::ExtMonUCITof).c_str();
    }

    static char const * StoppingTarget(){
      return StepInstanceName::name(StepInstanceName::stoppingtarget).c_str();
    }

    static char const * ProductionTargetCoreSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetCoreSection).c_str();
    }
    static char const * ProductionTargetStartingCoreSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetStartingCoreSection).c_str();
    }
    static char const * ProductionTargetFinStartingSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetFinStartingSection).c_str();
    }

    static char const * ProductionTargetNegativeEndRing(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetNegativeEndRing).c_str();
    }

    static char const * ProductionTargetPositiveEndRing(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetPositiveEndRing).c_str();
    }
    static char const * ProductionTargetFinSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetFinSection).c_str();
    }
    static char const * ProductionTargetFinTopSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetFinTopSection).c_str();
    }
    static char const * ProductionTargetFinTopStartingSection(){
      return StepInstanceName::name(StepInstanceName::ProductionTargetFinTopStartingSection).c_str();
    }

    static char const * CRSScintillatorBar(){
      return StepInstanceName::name(StepInstanceName::CRV).c_str();
    }

    static char const * TrackerPlaneSupport(){
      return StepInstanceName::name(StepInstanceName::trackerDS).c_str();
    }

    static char const * ProtonAbsorber() {
      return StepInstanceName::name(StepInstanceName::protonabsorber).c_str();
    }

    static char const * PSVacuum() {
      return StepInstanceName::name(StepInstanceName::PSVacuum).c_str();
    }

    static char const * TrackerSWires(){
      return StepInstanceName::name(StepInstanceName::trackerSWires).c_str();
    }

    static char const * ITrackerFWires(){
      return StepInstanceName::name(StepInstanceName::itrackerFWires).c_str();
    }

    static char const * TrackerWalls(){
      return StepInstanceName::name(StepInstanceName::trackerWalls).c_str();
    }

    static char const * STMDet() {
      return StepInstanceName::name(StepInstanceName::STMDet).c_str();
    }

    static char const * panelEBKey() {
      return StepInstanceName::name(StepInstanceName::panelEBKey).c_str();
    }

    static char const * DSCableRun() {
      return StepInstanceName::name(StepInstanceName::DSCableRun).c_str();
    }

  };

} // namespace mu2e

#endif /* Mu2eG4_SensitiveDetectorName_hh */
