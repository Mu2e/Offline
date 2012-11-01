#ifndef ConditionsService_ExtMonFNALConditions_hh
#define ConditionsService_ExtMonFNALConditions_hh
//
// $Id: ExtMonFNALConditions.hh,v 1.2 2012/11/01 23:39:26 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:39:26 $
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

    double t0() const { return t0_; }

    // absolute (K)
    double temperature() const { return temperature_; }

    double biasVoltage() const { return biasVoltage_; }

  private:
    int numClockTicksPerDebuncherPeriod_;
    double clockTick_;
    double t0_;

    double temperature_;
    double biasVoltage_;
  };

}

#endif /* ConditionsService_ExtMonFNALConditions_hh */
