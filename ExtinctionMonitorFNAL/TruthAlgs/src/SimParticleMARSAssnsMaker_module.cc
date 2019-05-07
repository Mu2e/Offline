// Associate SimParticles to MARSInfo in a newly produced compressed MARSInfoCollection.
//
// $Id: SimParticleMARSAssnsMaker_module.cc,v 1.3 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author Andrei Gaponenko
//

#include <string>
#include <memory>
#include <iostream>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleMARSAssns.hh"

#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class SimParticleMARSAssnsMaker : public art::EDProducer {

    public:
      explicit SimParticleMARSAssnsMaker(fhicl::ParameterSet const& pset)
        : EDProducer{pset}
        , simParticlesModuleLabel_(pset.get<std::string>("simParticlesModuleLabel"))
        , simParticlesInstanceName_(pset.get<std::string>("simParticlesInstanceName", ""))
        , marsInfoModuleLabel_(pset.get<std::string>("marsInfoModuleLabel"))
        , marsInfoInstanceName_(pset.get<std::string>("marsInfoInstanceName", ""))
      {
        produces<MARSInfoCollection>();
        produces<SimParticleMARSAssns>();
      }

      virtual void produce(art::Event& evt);

    private:
      std::string simParticlesModuleLabel_;
      std::string simParticlesInstanceName_;
      std::string marsInfoModuleLabel_;
      std::string marsInfoInstanceName_;
    };

    //================================================================
    void SimParticleMARSAssnsMaker::produce(art::Event& event) {

      std::unique_ptr<MARSInfoCollection> info(new MARSInfoCollection());
      std::unique_ptr<SimParticleMARSAssns> assns(new SimParticleMARSAssns());

      const art::ProductID infoPID = event.getProductID<MARSInfoCollection>();
      const art::EDProductGetter *infoGetter = event.productGetter(infoPID);

      art::Handle<GenParticleCollection> genpartsh;
      event.getByLabel(marsInfoModuleLabel_, marsInfoInstanceName_, genpartsh);
      art::FindOne<MARSInfo> infoFinder(genpartsh, event, art::InputTag(marsInfoModuleLabel_, marsInfoInstanceName_));

      art::Handle<SimParticleCollection> inpartsh;
      event.getByLabel(simParticlesModuleLabel_, simParticlesInstanceName_, inpartsh);
      const SimParticleCollection& inparts(*inpartsh);

      // Compute all primaries corresponding to the input collection
      typedef std::set<art::Ptr<SimParticle> > PrimarySet;
      PrimarySet primaries;
      SimParticleParentGetter pg(event);
      for(SimParticleCollection::const_iterator i=inparts.begin(); i!=inparts.end(); ++i) {

        art::Ptr<SimParticle> current(inpartsh, i->first.asInt());
        art::Ptr<SimParticle> next(pg.parent(current));

        while(next) {
          current = next;
          next = pg.parent(current);
        }

        primaries.insert(current);
      }

      // Record infos for each primary
      for(PrimarySet::const_iterator i=primaries.begin(); i!=primaries.end(); ++i) {

        art::Ptr<GenParticle> gp = (*i)->genParticle();
        if(!gp) {
          throw cet::exception("BADINPUTS")<<"Primary SimParticle has does not have genParticle\n";
        }

        const MARSInfo& mi = infoFinder.at(gp.key()).ref();
        info->push_back(mi);

        assns->addSingle(*i, art::Ptr<MARSInfo>(infoPID, info->size()-1, infoGetter));
      }

      //----------------
      event.put(std::move(info));
      event.put(std::move(assns));
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::SimParticleMARSAssnsMaker)
