#include "Offline/ConditionsService/inc/ExtMonFNALConditions.hh"

#include "cetlib_except/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

namespace mu2e {

  ExtMonFNALConditions::ExtMonFNALConditions(const SimpleConfig& config)
    : numClockTicksPerDebuncherPeriod_(config.getInt("extMonFNAL.numClockTicksPerDebuncherPeriod"))
    , clockTick_(GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod()/numClockTicksPerDebuncherPeriod_)
    , temperature_(config.getDouble("extMonFNAL.temperature")*CLHEP::kelvin)
    , biasVoltage_(config.getDouble("extMonFNAL.biasVoltage")*CLHEP::volt)
  {
    if(config.hasName("extMonFNAL.clockTick")) {
      throw cet::exception("CONFIG")<<"extMonFNAL.clockTick parameter should not be used.  Superseded by numClockTicksPerDebuncherPeriod.\n";
    }
  }

}
