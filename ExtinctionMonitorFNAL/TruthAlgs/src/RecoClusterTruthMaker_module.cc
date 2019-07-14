// Associate truth info to ExtMonFNALRecoClusters.
//
// $Id: RecoClusterTruthMaker_module.cc,v 1.3 2013/03/15 15:52:04 kutschke Exp $
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

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"

#include "MCDataProducts/inc/SimParticle.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      typedef std::map<art::Ptr<SimParticle>, double> ClusterCharges;
    }

    //================================================================
    class RecoClusterTruthMaker : public art::EDProducer {

    public:
      explicit RecoClusterTruthMaker(fhicl::ParameterSet const& pset)
        : EDProducer{pset}
        , recoClusterModuleLabel_(pset.get<std::string>("recoClusterModuleLabel"))
        , recoClusterInstanceName_(pset.get<std::string>("recoClusterInstanceName", ""))
        , hitTruthModuleLabel_(pset.get<std::string>("hitTruthModuleLabel"))
        , hitTruthInstanceName_(pset.get<std::string>("hitTruthInstanceName", ""))
      {
        produces<ExtMonFNALRecoClusterTruthAssn>();
      }

      virtual void produce(art::Event& evt);

    private:
      std::string recoClusterModuleLabel_;
      std::string recoClusterInstanceName_;
      std::string hitTruthModuleLabel_;
      std::string hitTruthInstanceName_;

      void getCharges(ClusterCharges *out, const art::Event& event, const ExtMonFNALRecoCluster& c);
    };

    //================================================================
    void RecoClusterTruthMaker::produce(art::Event& event) {

      std::unique_ptr<ExtMonFNALRecoClusterTruthAssn> outTruth(new ExtMonFNALRecoClusterTruthAssn());

      art::Handle<ExtMonFNALRecoClusterCollection> hrecoClusters;
      event.getByLabel(recoClusterModuleLabel_, recoClusterInstanceName_, hrecoClusters);
      const ExtMonFNALRecoClusterCollection& recoClusters(*hrecoClusters);

      for(unsigned ir = 0; ir < recoClusters.size(); ++ir) {
        ClusterCharges cc;
        getCharges(&cc, event, recoClusters[ir]);
        for(ClusterCharges::const_iterator i = cc.begin(); i != cc.end(); ++i) {
          outTruth->addSingle(i->first,
                              art::Ptr<ExtMonFNALRecoCluster>(hrecoClusters, ir),
                              ExtMonFNALRecoClusterTruthBits(i->second)
                              );
        }
      }

      event.put(std::move(outTruth));
    }

    //================================================================
    void RecoClusterTruthMaker::getCharges(ClusterCharges *out, const art::Event& event, const ExtMonFNALRecoCluster& c) {

      typedef ExtMonFNALRawCluster::Hits Hits;
      const Hits& rawHits = c.raw()->hits();

      art::FindManyP<SimParticle,ExtMonFNALHitTruthBits>
        fp(rawHits, event, art::InputTag(hitTruthModuleLabel_, hitTruthInstanceName_));

      for(unsigned ihit = 0; ihit < rawHits.size(); ++ihit) {

        std::vector<art::Ptr<SimParticle> > particles;
        std::vector<const ExtMonFNALHitTruthBits*> charges;

        const unsigned np = fp.get(ihit, particles, charges);
        for(unsigned ip=0; ip<np; ++ip) {
          (*out)[particles[ip]] += charges[ip]->charge();
        }
      }
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::RecoClusterTruthMaker)
