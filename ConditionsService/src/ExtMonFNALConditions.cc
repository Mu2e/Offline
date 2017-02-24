#include "ConditionsService/inc/ExtMonFNALConditions.hh"

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "ConditionsService/inc/AcceleratorParams.hh"

namespace mu2e {

  ExtMonFNALConditions::ExtMonFNALConditions(const AcceleratorParams& accp,
                                             const SimpleConfig& config)
    : numClockTicksPerDebuncherPeriod_(config.getInt("extMonFNAL.numClockTicksPerDebuncherPeriod"))
    , clockTick_(accp.deBuncherPeriod/numClockTicksPerDebuncherPeriod_)
    , temperature_(config.getDouble("extMonFNAL.temperature")*CLHEP::kelvin)
    , biasVoltage_(config.getDouble("extMonFNAL.biasVoltage")*CLHEP::volt)
  {
    if(config.hasName("extMonFNAL.clockTick")) {
      throw cet::exception("CONFIG")<<"extMonFNAL.clockTick parameter should not be used.  Superseded by numClockTicksPerDebuncherPeriod.\n";
    }
  }

}
