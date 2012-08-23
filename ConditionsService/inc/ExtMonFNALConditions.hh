#ifndef ConditionsService_ExtMonFNALConditions_hh
#define ConditionsService_ExtMonFNALConditions_hh
//
// $Id: ExtMonFNALConditions.hh,v 1.1 2012/08/23 23:41:52 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:41:52 $
//
// Original author Andrei Gaponenko
//

#include <iostream>

#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

namespace mu2e
{
  class SimpleConfig;

  struct ExtMonFNALConditions: virtual public ConditionsEntity {

    explicit ExtMonFNALConditions(const SimpleConfig& config);

    // Time bin width, ns
    double clockTick() const { return clockTick_; }

    double t0() const { return t0_; }

    // absolute (K)
    double temperature() const { return temperature_; }

    double biasVoltage() const { return biasVoltage_; }

  private:
    double clockTick_;
    double t0_;

    double temperature_;
    double biasVoltage_;
  };

}

#endif /* ConditionsService_ExtMonFNALConditions_hh */
