#ifndef ConditionsService_ExtMonFNALConditions_hh
#define ConditionsService_ExtMonFNALConditions_hh
//
//
// Original author Andrei Gaponenko
//

#include <iostream>

#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;
  class AcceleratorParams;

  struct ExtMonFNALConditions: virtual public ConditionsEntity {

    explicit ExtMonFNALConditions(const AcceleratorParams& accp, const SimpleConfig& config);

    int numClockTicksPerDebuncherPeriod() const { return numClockTicksPerDebuncherPeriod_; }

    // Time bin width, ns
    double clockTick() const { return clockTick_; }

    // absolute (K)
    double temperature() const { return temperature_; }

    double biasVoltage() const { return biasVoltage_; }

  private:
    int numClockTicksPerDebuncherPeriod_;
    double clockTick_;

    double temperature_;
    double biasVoltage_;
  };

}

#endif /* ConditionsService_ExtMonFNALConditions_hh */
