//  Two-way conversion between module id and sequential ("dense") module number.
//
// Evan Schiewe, 2013

#ifndef EXTMONFNALMODULEIDCONVERTER_HH
#define EXTMONFNALMODULEIDCONVERTER_HH
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "DataProducts/inc/ExtMonFNALModuleDenseId.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

namespace mu2e {

  class ExtMonFNALModuleIdConverter {
  public:

    // should never be called with invalid id
    ExtMonFNALModuleDenseId denseModuleNumber(const ExtMonFNALModuleId& id) const;

    // should never be called with invalid dense id
    ExtMonFNALModuleId moduleId(ExtMonFNALModuleDenseId did) const;

    ExtMonFNALModuleDenseId getModuleDenseId(unsigned plane, unsigned mod) const;
    ExtMonFNALModuleId getModuleId(unsigned plane, unsigned mod) const;

    explicit ExtMonFNALModuleIdConverter(const mu2e::ExtMonFNAL::ExtMon& extmon);

  private:
    const mu2e::ExtMonFNAL::ExtMon *extmon_;
  };

} // namespace mu2e

#endif/*EXTMONFNALMODULEIDCONVERTER_HH*/
