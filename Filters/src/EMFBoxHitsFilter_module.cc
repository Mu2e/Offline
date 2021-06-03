// Filters outputs of ExtMonFNAL "room" jobs to select VD hits and
// stopped muons that should be passed to the next stage "box" jobs.
//
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {

      // Use art::Ptr instead of bare ptr to get the desired sorting
      typedef std::set<art::Ptr<SimParticle> > TrackSet;

      // Adapter for compressSimParticleCollection()
      class ParticleSelector {
      public:
        ParticleSelector(const TrackSet& m) {
          for(TrackSet::const_iterator i = m.begin(); i!=m.end(); ++i) {
            m_keys.insert((*i)->id());
          }
        }

        bool operator[]( cet::map_vector_key key ) const {
          return m_keys.find(key) != m_keys.end();
        }

      private:
        std::set<cet::map_vector_key> m_keys;
      };

    }

    //================================================================
    class EMFBoxHitsFilter : public art::EDFilter {
      std::string hitsModuleLabel_;
      std::string hitsInstanceName_;
      std::string simParticlesModuleLabel_;
      std::string simParticlesInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      double cutEKineAtStop_;
      std::vector<double> cutStoppedMuonPosition_;

      // statistics
      unsigned nVDHits_;
      unsigned nStoppedMuons_;
      unsigned nPassedEvents_;

      const ExtMon *extmon_;

    public:
      explicit EMFBoxHitsFilter(const fhicl::ParameterSet& pset);
      virtual bool beginRun(art::Run& run);
      virtual bool filter(art::Event& event);
      virtual void endJob();
    };

    //================================================================
    EMFBoxHitsFilter::EMFBoxHitsFilter(const fhicl::ParameterSet& pset)
      : art::EDFilter{pset}
      , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
      , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))
      , simParticlesModuleLabel_(pset.get<std::string>("simParticlesModuleLabel"))
      , simParticlesInstanceName_(pset.get<std::string>("simParticlesInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
      , cutEKineAtStop_(pset.get<double>("cutEKineAtStop"))
      , cutStoppedMuonPosition_(pset.get<std::vector<double> >("cutStoppedMuonPosition"))


      , nVDHits_()
      , nStoppedMuons_()
      , nPassedEvents_()

      , extmon_()
    {
      produces<StepPointMCCollection>();
      produces<SimParticleCollection>();
    }

    //================================================================
    bool EMFBoxHitsFilter::beginRun(art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        extmon_ = &*extmon;
      }
      else {
        GeomHandle<ExtMon> extmon;
        extmon_ = &*extmon;
      }
      return true;
    }

    //================================================================
    bool EMFBoxHitsFilter::filter(art::Event& event) {

      std::unique_ptr<StepPointMCCollection> outhits(new StepPointMCCollection());
      std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());

      //----------------
      // First, find all stopped muons that we wish to save

      art::Handle<SimParticleCollection> inparticlesh;
      event.getByLabel(simParticlesModuleLabel_, simParticlesInstanceName_, inparticlesh);
      const SimParticleCollection& inparticles(*inparticlesh);

      TrackSet stoppedMuons;

      for(SimParticleCollection::const_iterator i=inparticles.begin(); i!=inparticles.end(); ++i) {
        const SimParticle& sp = i->second;
        if(std::abs(sp.pdgId())==13) {
          const double eKine = sp.endMomentum().e() - sp.endMomentum().m();
          if(eKine < cutEKineAtStop_) {

            const CLHEP::Hep3Vector extMonPos = extmon_->mu2eToExtMon_position(sp.endPosition());

            // Only record stops in the vicinity of ExtMonFNAL
            if((std::abs(extMonPos.x()) < cutStoppedMuonPosition_[0]) &&
               (std::abs(extMonPos.y()) < cutStoppedMuonPosition_[1]) &&
               (std::abs(extMonPos.z()) < cutStoppedMuonPosition_[2])
               ) {
              stoppedMuons.insert(art::Ptr<SimParticle>(inparticlesh, i->first.asUint()));
              ++nStoppedMuons_;
            }

          } // if(stopped)
        } // if(muon)
      } // for(particles)

      //----------------
      // Now select VD hits we want, skipping hits from stopped muons
      // that would be double counted

      art::Handle<StepPointMCCollection> inhitsh;
      event.getByLabel(hitsModuleLabel_, hitsInstanceName_, inhitsh);
      const StepPointMCCollection& inhits(*inhitsh);

      TrackSet particlesWithHits;
      for(StepPointMCCollection::const_iterator i=inhits.begin(); i!=inhits.end(); ++i) {
        switch(i->volumeId()) {
        case VirtualDetectorId::EMFBoxFront: case VirtualDetectorId::EMFBoxBack:
        case VirtualDetectorId::EMFBoxSW: case VirtualDetectorId::EMFBoxNE:
        case VirtualDetectorId::EMFBoxBottom: case VirtualDetectorId::EMFBoxTop:

          if(stoppedMuons.find(i->simParticle()) == stoppedMuons.end()) { // not a stopped muon
            outhits->push_back(*i);
            particlesWithHits.insert(i->simParticle());
            ++nVDHits_;
          }

          break;

        default:
          break;
        }
      }

      //----------------
      // Prepare a complete set of particles to keep, including all the mothers
      particlesWithHits.insert(stoppedMuons.begin(), stoppedMuons.end());
      SimParticleParentGetter pg(event);
      TrackSet particlesToKeep;
      for(TrackSet::const_iterator i=particlesWithHits.begin(); i!=particlesWithHits.end(); ++i) {
        art::Ptr<SimParticle> particle(*i);
        // Get all the parents
        do {
          particlesToKeep.insert(particle);
          particle = pg.parent(particle);
        } while(particle);
      }

      //----------------
      // Prepare the output particle collection
      ParticleSelector selector(particlesToKeep);
      art::ProductID newParticlesPID(event.getProductID<SimParticleCollection>());
      const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));
      compressSimParticleCollection(newParticlesPID, newParticlesGetter, *inparticlesh, selector, *outparts);

      //----------------
      // Updates pointes in the hit collection to point to the new SimParticleCollection
      for(StepPointMCCollection::iterator i=outhits->begin(); i!=outhits->end(); ++i) {
        art::Ptr<SimParticle> oldPtr(i->simParticle());
        i->simParticle() = art::Ptr<SimParticle>(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
      }

      //----------------
      const bool passed = !outparts->empty();
      if(passed) {
        ++nPassedEvents_;
      }
      event.put(std::move(outparts));
      event.put(std::move(outhits));
      return passed;;
    }

    //================================================================
    void EMFBoxHitsFilter::endJob() {
      mf::LogInfo("Summary")
        <<"EMFBoxHitsFilter_module: Number of events passing the filter: "<<nPassedEvents_
        <<", nVDHits = "<<nVDHits_
        <<", nStoppedMuons = "<<nStoppedMuons_
        << "\n";
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFBoxHitsFilter);
