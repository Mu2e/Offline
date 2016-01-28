//
// The btrk code needs some information and tools that are supplied externally.
// This service is the agent that supplies the that information and those tools.
//
// Also supply the reconstruction geometry used by BTrk.
//
// Original author Rob Kutschke
//
// Notes:
// 1) For correct results, the beginRun member functions of GeometryService and
//    ConditionsService must be called prior to calling the beginRun member function
//    of this service.  To force this to happen, get service handles to those
//    services in the c'tor to this service and do so prior to registering callbacks
//    for this service.

#include "BTrkHelper/inc/BTrkHelper.hh"

#include "ConditionsService/inc/ConditionsService.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "BTrk/BaBar/ExternalInfo.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"

mu2e::BTrkHelper::BTrkHelper(fhicl::ParameterSet const& pset,
                             art::ActivityRegistry&     registry):
  particleInfo_(),
  fileFinder_(pset.get<fhicl::ParameterSet>("dictionaries")),
  detectorModelPSet_(pset.get<fhicl::ParameterSet>("Mu2eDetectorModel"))
{

  // See note 1).
  art::ServiceHandle<GeometryService> geom [[gnu::unused]];
  art::ServiceHandle<ConditionsService> conditions [[gnu::unused]];

  // Register callbacks.
  registry.sPreBeginRun.watch   ( this, &BTrkHelper::beginRun    );
  registry.sPreBeginSubRun.watch( this, &BTrkHelper::beginSubRun );

  // These are available at c'tor time - so push them to the BaBar code now.
  ExternalInfo::set( &fileFinder_);
  ExternalInfo::set( &particleInfo_);
}

void
mu2e::BTrkHelper::beginRun(art::Run const &) {

  // For now we assume that the geometry does not change over the course of a job.
  // FIXME: later on allow geometry to change.
  // FIXME: do we want to remove compile time dependence on TTracker and Mu2eDetectorModel
  //        by holding Mu2eDetectorModel as pointer to a base class and by using factory function?
  if ( detectorModel_.get() == nullptr ){
    GeomHandle<TTracker> tt;
    detectorModel_ = std::make_unique<Mu2eDetectorModel>( detectorModelPSet_, *tt );
  }
}

void
mu2e::BTrkHelper::beginSubRun(art::SubRun const &) {
  // When we have information that needs to be pushed at beginSubRun time, then this is where to do it.
}

DEFINE_ART_SERVICE(mu2e::BTrkHelper);
