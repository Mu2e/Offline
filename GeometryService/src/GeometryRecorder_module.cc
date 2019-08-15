// Attach select GeometryService data to the run object.
//
// Andrei Gaponenko, 2012

#include <memory>

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"
#include "Mu2eHallGeom/inc/Mu2eHall.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

namespace mu2e {

  //================================================================
  class GeometryRecorder : public art::EDProducer {
    template<class DET> void put(art::Run& run) {
      run.put(std::unique_ptr<DET>(new DET(*GeomHandle<DET>())));
    }

  public:
    explicit GeometryRecorder(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event&) {/*No per-event actions*/}
    virtual void endRun(art::Run&);
  };

  //================================================================
  GeometryRecorder::GeometryRecorder(const fhicl::ParameterSet& pset) :
    art::EDProducer{pset} 
 {
    produces<Mu2eHall, art::InRun>();
    produces<ProductionTarget, art::InRun>();
    produces<ProtonBeamDump, art::InRun>();
    produces<ExtMonFNALBuilding, art::InRun>();
    produces<ExtMonFNAL::ExtMon, art::InRun>();
    produces<PSEnclosure, art::InRun>();
    produces<PSShield, art::InRun>();
    produces<PSVacuum, art::InRun>();
    produces<Mu2eEnvelope, art::InRun>();

    produces<WorldG4, art::InRun>();
  }

  //================================================================
  // Better do the work in endRun() than beginRun(), otherwise
  // the WorldG4 object may not yet be available.
  void GeometryRecorder::endRun(art::Run& run) {
    put<Mu2eHall>(run);
    put<ProductionTarget>(run);
    put<ProtonBeamDump>(run);
    put<ExtMonFNALBuilding>(run);
    put<ExtMonFNAL::ExtMon>(run);
    put<PSEnclosure>(run);
    put<PSShield>(run);
    put<PSVacuum>(run);
    put<Mu2eEnvelope>(run);

    art::ServiceHandle<GeometryService> geom;
    if(geom->hasElement<WorldG4>()) {
      put<WorldG4>(run);
    }
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GeometryRecorder);
