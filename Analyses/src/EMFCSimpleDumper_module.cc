// Ntuple dumper for ExtMon.  The "simple" version writes hits for a single VD.
//
// Andrei Gaponenko, 2012, 2015

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// ROOT
#include "TDirectory.h"
#include "TH1.h"
#include "TTree.h"

// externals
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

// std
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <limits>


#include <iostream>
//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff) do {} while(0)

namespace mu2e {
  namespace ExtMonFNAL {

    class Mu2eCoordinates {};

    //================================================================
    double getCharge(PDGCode::type pdgId) {
      // unlike generic conditions, MC particle data
      // should not change run-to-run, so static is safe
      // use static for efficiency
      static GlobalConstantsHandle<ParticleDataTable> pdt;

      ParticleDataTable::maybe_ref info = pdt->particle(pdgId);

      // Negatives will be rejected by the analysis
      // Use an impossible number which is in the range of the "charge" histogram to signal "uninitialized"
      double charge(-0.2);

      if(!info.isValid()) {
        std::cout<<"AG: warning: no valid PDG info for pdgId = "<<pdgId<<" using charge = "<<charge<<std::endl;
      }
      else {
        charge = info.ref().charge();
      }

      return charge;
    }

    //================================================================
    struct VDHit {
      float pmag;
      float px;
      float py;
      float pz;

      float x;
      float y;
      float z;

      float time;

      float charge;
      int   pdgId;
      unsigned particleId;

      VDHit() : pmag(std::numeric_limits<double>::quiet_NaN())
              , px(std::numeric_limits<double>::quiet_NaN())
              , py(std::numeric_limits<double>::quiet_NaN())
              , pz(std::numeric_limits<double>::quiet_NaN())

              , x(std::numeric_limits<double>::quiet_NaN())
              , y(std::numeric_limits<double>::quiet_NaN())
              , z(std::numeric_limits<double>::quiet_NaN())

              , time(std::numeric_limits<double>::quiet_NaN())

              , charge(std::numeric_limits<double>::quiet_NaN())
              , pdgId(0)
              , particleId(-1U)
      {}

      //----------------------------------------------------------------
      // Constructor using mu2e coords
      VDHit(const Mu2eCoordinates& , const StepPointMC& hit)
        : pmag(hit.momentum().mag())
        , px(hit.momentum().x())
        , py(hit.momentum().y())
        , pz(hit.momentum().z())

        , x(hit.position().x())
        , y(hit.position().y())
        , z(hit.position().z())

        , time(hit.time())
        , charge(getCharge(hit.simParticle()->pdgId()))

        , pdgId(hit.simParticle()->pdgId())
        , particleId(hit.simParticle()->id().asUint())
      {}

      //----------------------------------------------------------------
      // Constructor converting to the beam dump system
      VDHit(const ProtonBeamDump& dump, const StepPointMC& hit)
        : pmag(hit.momentum().mag())
        , px(), py(), pz(), x(), y(), z()
        , time(hit.time())
        , charge(getCharge(hit.simParticle()->pdgId()))
        , pdgId(hit.simParticle()->pdgId())
        , particleId(hit.simParticle()->id().asUint())
      {
        const CLHEP::Hep3Vector p = dump.mu2eToBeamDump_momentum(hit.momentum());
        px = p.x();
        py = p.y();
        pz = p.z();

        const CLHEP::Hep3Vector v = dump.mu2eToBeamDump_position(hit.position());
        x = v.x();
        y = v.y();
        z = v.z();
      }

      //----------------------------------------------------------------
      // Constructor converting to the ExtMonFNAL detector system
      VDHit(const ExtMonFNALPlaneStack& det, const StepPointMC& hit)
        : pmag(hit.momentum().mag())
        , px(), py(), pz(), x(), y(), z()
        , time(hit.time())
        , charge(getCharge(hit.simParticle()->pdgId()))
        , pdgId(hit.simParticle()->pdgId())
        , particleId(hit.simParticle()->id().asUint())
      {
        const CLHEP::Hep3Vector p = det.mu2eToStack_momentum(hit.momentum());
        px = p.x();
        py = p.y();
        pz = p.z();

        const CLHEP::Hep3Vector v = det.mu2eToStack_position(hit.position());
        x = v.x();
        y = v.y();
        z = v.z();
      }

    }; // struct VDHit

    //================================================================
    class EMFCSimpleDumper : public art::EDAnalyzer {
      std::string triggerInputCollection_;
      unsigned triggerVD_;

      // Coordinate conversion
      const ProtonBeamDump *pdump_;
      const ExtMon *pdet_;

      // Members needed to write the ntuple
      TTree *nt_;

      VDHit trigHit_;

      // a helper to decide what coordinate conversion to use
      VDHit makeVDHit(const StepPointMC& hit) const;

    public:
      explicit EMFCSimpleDumper(const fhicl::ParameterSet& pset);
      virtual void beginJob();
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFCSimpleDumper::EMFCSimpleDumper(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , triggerInputCollection_(pset.get<std::string>("triggerInputCollection"))
      , triggerVD_(pset.get<unsigned>("triggerVD"))
      , pdump_(0)
      , pdet_(0)
      , nt_(0)
    {}

    //================================================================
    void EMFCSimpleDumper::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;
      static const char branchDesc[] = "p/F:px/F:py/F:pz/F:x/F:y/F:z/F:time/F:charge/F:pdgId/I:particleId/i";
      nt_ = tfs->make<TTree>( "nt", "EMFCSimpeDumper ntuple");
      nt_->Branch("trig", &trigHit_, branchDesc);
    }

    //================================================================
    void EMFCSimpleDumper::beginRun(const art::Run& run) {
      // // Note: it may be useful to retrieve geometry objects from the
      // // input file, as opposed to using GeometryService that might have
      // // a different version than used to produce our input data.
      //
      // art::Handle<ProtonBeamDump> dumph;
      // run.getByLabel(geomModuleLabel_, "", dumph);
      // pdump_ = &*dumph;
      pdump_ = GeomHandle<ProtonBeamDump>().get();

      std::cout<<"EMFCSimpleDumper::beginRun(): got ProtonBeamDump: coreCenterInMu2e = "
               <<pdump_->coreCenterInMu2e()
               <<", coreRotY = "<<pdump_->coreRotY()/CLHEP::degree<<" degrees"
               <<std::endl;

      // art::Handle<ExtMon> deth;
      // run.getByLabel(geomModuleLabel_, "", deth);
      // pdet_ = &*deth;
      pdet_ = GeomHandle<ExtMon>().get();

      std::cout<<"EMFCSimpleDumper::beginRun(): got ExtMonFNAL: detectorCenterInMu2e = "
               <<pdet_->detectorCenterInMu2e()
               <<", rot = "<<pdet_->detectorRotationInMu2e()
               <<std::endl;
    }

    //================================================================
    VDHit EMFCSimpleDumper::makeVDHit(const StepPointMC& hit) const {
      switch(hit.virtualDetectorId()) {
      default:
        return VDHit(Mu2eCoordinates(), hit);
      case VirtualDetectorId::EMFC1Entrance: case VirtualDetectorId::EMFC1Exit:
        return VDHit(*pdump_, hit);
      case VirtualDetectorId::EMFC2Entrance: case VirtualDetectorId::EMFC2Exit:
      case VirtualDetectorId::EMFDetectorUpEntrance: case VirtualDetectorId::EMFDetectorUpExit:
        return VDHit(pdet_->up(), hit);
      case VirtualDetectorId::EMFDetectorDnEntrance: case VirtualDetectorId::EMFDetectorDnExit:
        return VDHit(pdet_->dn(), hit);
      }
    }

    //================================================================
    void EMFCSimpleDumper::analyze(const art::Event& event) {
      const auto& trighits = event.getValidHandle<StepPointMCCollection>(triggerInputCollection_);
      for(const auto& i: *trighits) {
        if(i.volumeId() == triggerVD_) {
          trigHit_ = makeVDHit(i);
          nt_->Fill();
        }
      }
    }

    //================================================================

  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFCSimpleDumper);
