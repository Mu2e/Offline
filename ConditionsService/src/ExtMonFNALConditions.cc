#include "ConditionsService/inc/ExtMonFNALConditions.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"


namespace mu2e {

  ExtMonFNALConditions::ExtMonFNALConditions(const SimpleConfig& config)
    : clockTick_(config.getDouble("extMonFNAL.clockTick")*CLHEP::ns)
    , t0_(config.getDouble("extMonFNAL.t0")*CLHEP::ns)
    , temperature_(config.getDouble("extMonFNAL.temperature")*CLHEP::kelvin)
    , biasVoltage_(config.getDouble("extMonFNAL.biasVoltage")*CLHEP::volt)
  {}

}
