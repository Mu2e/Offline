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

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

namespace mu2e {
  namespace ExtMonFNAL {

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

    //================================================================
    class EMFPixelSimFilter : public art::EDFilter {
      std::string hitsModuleLabel_;
      std::string simParticlesModuleLabel_;

      unsigned cutMinPlanes_;

      unsigned nPassed_;

    public:
      explicit EMFPixelSimFilter(const fhicl::ParameterSet& pset);
      virtual bool filter(art::Event& event);
      virtual void endJob();
    };

    //================================================================
    EMFPixelSimFilter::EMFPixelSimFilter(const fhicl::ParameterSet& pset)
      : art::EDFilter{pset}
      , hitsModuleLabel_(pset.get<std::string>("simHitsModuleLabel"))
      , simParticlesModuleLabel_(pset.get<std::string>("simParticlesModuleLabel"))
      , cutMinPlanes_(pset.get<unsigned>("cutMinPlanes"))
      , nPassed_()
    {
      produces<SimParticleCollection>();
      produces<ExtMonFNALSimHitCollection>();
    }

    //================================================================
    bool EMFPixelSimFilter::filter(art::Event& event) {
      std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());
      std::unique_ptr<ExtMonFNALSimHitCollection> outhits(new ExtMonFNALSimHitCollection());

      std::set<unsigned> hitPlanes;
      SimParticleParentGetter pg(event);

      art::Handle<ExtMonFNALSimHitCollection> inhitsh;
      event.getByLabel(hitsModuleLabel_, inhitsh);
      const ExtMonFNALSimHitCollection& inhits(*inhitsh);

      // Use art::Ptr instead of bare ptr to get the desired sorting
      TrackSet particlesWithHits;
      for(ExtMonFNALSimHitCollection::const_iterator i=inhits.begin(); i != inhits.end(); ++i) {
        // Record sim particle making the hit, and all its parents
        art::Ptr<SimParticle> particle = i->simParticle();
        do {
          particlesWithHits.insert(particle);
          particle = pg.parent(particle);
        } while(particle);

        // Record hit plane
        hitPlanes.insert(i->moduleId().plane());
      }

      const bool passed = (cutMinPlanes_ <= hitPlanes.size());
      if(passed) {
        ++nPassed_;

        art::Handle<SimParticleCollection> inparts;
        event.getByLabel(simParticlesModuleLabel_, inparts);

        ParticleSelector selector(particlesWithHits);

        art::ProductID newParticlesPID(event.getProductID<SimParticleCollection>());
        const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

        compressSimParticleCollection(newParticlesPID, newParticlesGetter, *inparts, selector, *outparts);
        //----------------
        // Copy all hits to the output, while updating their SimParticle pointers
        outhits->reserve(inhits.size());
        for(ExtMonFNALSimHitCollection::const_iterator i=inhits.begin(); i != inhits.end(); ++i) {
          art::Ptr<SimParticle> oldPtr(i->simParticle());
          art::Ptr<SimParticle> newPtr(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
          outhits->push_back(ExtMonFNALSimHit(i->moduleId(),
                                              newPtr,
                                              i->totalEnergyDeposit(),
                                              i->nonIonizingEnergyDeposit(),
                                              i->localStartPosition(),
                                              i->startTime(),
                                              i->localEndPosition(),
                                              i->endTime()));
        }
      } // if(passed)

      event.put(std::move(outparts));
      event.put(std::move(outhits));
      return passed;
    }

    //================================================================
    void EMFPixelSimFilter::endJob() {
      mf::LogInfo("Summary")
        << "EMFPixelSimFilter_module: Number of events passing the filter: "
        << nPassed_
        << "\n";
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFPixelSimFilter);
