// Evan Schiewe, 2013

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

#include <cassert>
namespace mu2e {

  //================================================================
  ExtMonFNALModuleIdConverter::ExtMonFNALModuleIdConverter(const mu2e::ExtMonFNAL::ExtMon& extmon)
    : extmon_(&extmon)
   {}

  //================================================================
  ExtMonFNALModuleDenseId ExtMonFNALModuleIdConverter::denseModuleNumber(const ExtMonFNALModuleId& id) const {
    return ExtMonFNALModuleDenseId(getModuleDenseId(id.plane(),id.number()));
  }

  //================================================================
  ExtMonFNALModuleId ExtMonFNALModuleIdConverter::moduleId(ExtMonFNALModuleDenseId did) const {
    unsigned int module = did.number();
    unsigned int plane = 0;

    while ( module >= extmon_->plane(plane).nModules() ) {
      module -= extmon_->plane(plane).nModules();
      plane++;
    }

    return ExtMonFNALModuleId(plane, module);
  }


 

  ExtMonFNALModuleDenseId ExtMonFNALModuleIdConverter::getModuleDenseId(unsigned plane, unsigned mod) const{
    unsigned res = mod;
    for (unsigned iplane = 0; iplane < plane; ++iplane)
      res += extmon_->plane(iplane).nModules();
    return ExtMonFNALModuleDenseId(res);
  }
  
  ExtMonFNALModuleId ExtMonFNALModuleIdConverter::getModuleId(unsigned plane, unsigned mod) const{
    ExtMonFNALModuleDenseId did = this->getModuleDenseId(plane, mod);
    return this->moduleId(did);
  }
  
} // namespace mu2e
