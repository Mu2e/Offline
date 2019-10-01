// Associate truth to track finding output.
//
// $Id: PatRecTruthMaker_module.cc,v 1.5 2013/03/15 15:52:04 kutschke Exp $
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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      struct TrackParticles {
        std::vector<art::Ptr<SimParticle> > particles;
        std::vector<unsigned> nCommonClusters;
      };
    }

    //================================================================
    class PatRecTruthMaker : public art::EDProducer {

    public:
      explicit PatRecTruthMaker(fhicl::ParameterSet const& pset)
        : EDProducer{pset}
        , patRecModuleLabel_(pset.get<std::string>("patRecModuleLabel"))
        , patRecInstanceName_(pset.get<std::string>("patRecInstanceName", ""))
        , clusterTruthModuleLabel_(pset.get<std::string>("clusterTruthModuleLabel"))
        , clusterTruthInstanceName_(pset.get<std::string>("clusterTruthInstanceName", ""))
      {
        produces<ExtMonFNALPatRecTruthAssns>();
      }

      virtual void produce(art::Event& evt);

    private:
      std::string patRecModuleLabel_;
      std::string patRecInstanceName_;
      std::string clusterTruthModuleLabel_;
      std::string clusterTruthInstanceName_;

      void getParticles(TrackParticles *out,
                        const art::Event& event,
                        const ExtMonFNALTrkFit& track);
    };

    //================================================================
    void PatRecTruthMaker::produce(art::Event& event) {

      std::unique_ptr<ExtMonFNALPatRecTruthAssns> outTruth(new ExtMonFNALPatRecTruthAssns());

      art::Handle<ExtMonFNALTrkFitCollection> htracks;
      event.getByLabel(patRecModuleLabel_, patRecInstanceName_, htracks);
      const ExtMonFNALTrkFitCollection& tracks(*htracks);

      for(unsigned itrack = 0; itrack < tracks.size(); ++itrack) {

        TrackParticles cc;
        getParticles(&cc, event, tracks[itrack]);

        art::FindManyP<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>
          particleClusterFinder(cc.particles,
                                event,
                                art::InputTag(clusterTruthModuleLabel_, clusterTruthInstanceName_));

        for(unsigned iparticle = 0; iparticle < cc.particles.size(); ++iparticle) {

          unsigned nParticleClusters = particleClusterFinder.at(iparticle).size();

          outTruth->addSingle(cc.particles[iparticle],
                              art::Ptr<ExtMonFNALTrkFit>(htracks, itrack),
                              ExtMonFNALTrkMatchInfo(cc.nCommonClusters[iparticle],
                                                     tracks[itrack].clusters().size(),
                                                     nParticleClusters)
                              );
        }
      }

      event.put(std::move(outTruth));
    }

    //================================================================
    void PatRecTruthMaker::getParticles(TrackParticles *out,
                                        const art::Event& event,
                                        const ExtMonFNALTrkFit& track)
    {
      art::FindManyP<SimParticle,ExtMonFNALRecoClusterTruthBits>
        simParticleFinder(track.clusters(),
                          event,
                          art::InputTag(clusterTruthModuleLabel_, clusterTruthInstanceName_));

      // particle to number of common hits
      typedef std::map<art::Ptr<SimParticle>, unsigned> TrackParticleMap;
      TrackParticleMap pm;

      for(unsigned icluster = 0; icluster < track.clusters().size(); ++icluster) {

        const std::vector<art::Ptr<SimParticle> >& particles = simParticleFinder.at(icluster);

        for(unsigned ip=0; ip<particles.size(); ++ip) {
          pm[particles[ip]] += 1;
        }

      }

      for(TrackParticleMap::const_iterator i = pm.begin(); i!=pm.end(); ++i) {
        out->particles.push_back(i->first);
        out->nCommonClusters.push_back(i->second);
      }
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::PatRecTruthMaker)
