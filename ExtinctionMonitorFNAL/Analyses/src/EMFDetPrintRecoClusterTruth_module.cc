// Printout ExtMonFNAL raw hits and associated truth
//
// Andrei Gaponenko, 2012

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/Ptr.h"

#include <iostream>
#include <string>
#include <vector>

namespace mu2e {

  //================================================================
  class EMFDetPrintRecoClusterTruth : public art::EDAnalyzer {
    std::string clusterModuleLabel_;  // emtpy label disables by-cluster printout
    std::string clusterInstanceName_;
    std::string truthModuleLabel_;
    std::string truthInstanceName_;
    std::string particleModuleLabel_;  // emtpy label disables by-particle printout
    std::string particleInstanceName_;

    unsigned cutMinClustersPerParticle_;

    void printByCluster(const art::Event& event);
    void printByParticle(const art::Event& event);

  public:
    explicit EMFDetPrintRecoClusterTruth(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetPrintRecoClusterTruth::EMFDetPrintRecoClusterTruth(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , clusterModuleLabel_(pset.get<std::string>("clusterModuleLabel"))
    , clusterInstanceName_(pset.get<std::string>("clusterInstanceName", ""))
    , truthModuleLabel_(pset.get<std::string>("truthModuleLabel"))
    , truthInstanceName_(pset.get<std::string>("truthInstanceName", ""))
    , particleModuleLabel_(pset.get<std::string>("particleModuleLabel"))
    , particleInstanceName_(pset.get<std::string>("particleInstanceName", ""))
    , cutMinClustersPerParticle_(pset.get<unsigned>("cutMinClustersPerParticle"))
  {}

  //================================================================
  void EMFDetPrintRecoClusterTruth::analyze(const art::Event& event) {
    if(!clusterModuleLabel_.empty()) {
      printByCluster(event);
    }
    if(!particleModuleLabel_.empty()) {
      printByParticle(event);
    }
  }

  //================================================================
  void EMFDetPrintRecoClusterTruth::printByCluster(const art::Event& event) {

    art::Handle<ExtMonFNALRecoClusterCollection> ih;
    event.getByLabel(clusterModuleLabel_, clusterInstanceName_, ih);

    const ExtMonFNALRecoClusterCollection& inputs(*ih);

    std::cout<<"EMFDetPrintRecoClusterTruth: by-cluster printout for event "<<event.id()<<std::endl;
    art::FindMany<SimParticle,ExtMonFNALRecoClusterTruthBits>
      r2t(ih, event, art::InputTag(truthModuleLabel_, truthInstanceName_));

    for(ExtMonFNALRecoClusterCollection::const_iterator i=inputs.begin(); i!=inputs.end(); ++i) {
      std::cout<<"    Cluster: "<<*i<<std::endl;

      const std::vector<const SimParticle*>& particles = r2t.at(i-inputs.begin());
      const std::vector<const ExtMonFNALRecoClusterTruthBits*>& charges = r2t.data(i-inputs.begin());

      std::cout<<"    got truth: num particles = "<<particles.size()<<", num charges "<<charges.size()<<std::endl;;
      for(unsigned ip=0; ip<particles.size(); ++ip) {
        std::cout<<"        particle id="<<particles[ip]->id()<<", charge "<<charges[ip]->charge()<<std::endl;
      }
    }
  }
  //================================================================
  void EMFDetPrintRecoClusterTruth::printByParticle(const art::Event& event) {

    std::cout<<"EMFDetPrintRecoClusterTruth: by-particle printout for event "<<event.id()
             <<", cutMinClustersPerParticle = "<<cutMinClustersPerParticle_
             <<std::endl;

    art::Handle<SimParticleCollection> ih;
    event.getByLabel(particleModuleLabel_, particleInstanceName_, ih);

    // FIXME: need to create sequence of SimParticles by hand
    // because map_vector does not work with FindMany
    // https://cdcvs.fnal.gov/redmine/issues/2967
    std::vector<art::Ptr<SimParticle> > particles;
    for(SimParticleCollection::const_iterator i = ih->begin(), iend = ih->end(); i != iend; ++i) {
      particles.push_back(art::Ptr<SimParticle>(ih, i->first.asUint()));
    }

    art::FindMany<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>
      t2r(particles, event, art::InputTag(truthModuleLabel_, truthInstanceName_));

    const SimParticleCollection& inputs(*ih);
    for(SimParticleCollection::const_iterator i=inputs.begin(); i!=inputs.end(); ++i) {

      const std::vector<const ExtMonFNALRecoCluster*>& clusters = t2r.at(i-inputs.begin());
      const std::vector<const ExtMonFNALRecoClusterTruthBits*>& charges = t2r.data(i-inputs.begin());

      if(clusters.size() >= cutMinClustersPerParticle_) {
        std::cout<<"    particle: "<<i->second.id()<<std::endl;
        std::cout<<"    got reco: num clusters = "<<clusters.size()<<", num charges "<<charges.size()<<std::endl;;
        for(unsigned ic=0; ic<clusters.size(); ++ic) {
          std::cout<<"        charge ="<<charges[ic]->charge()
                   <<", cluster = "<<*clusters[ic]<<", hits = { ";

          typedef ExtMonFNALRawCluster::Hits Hits;
          const Hits& raw = clusters[ic]->raw()->hits();
          for(Hits::const_iterator i = raw.begin(); i != raw.end(); ++i) {
            std::cout<<**i<<", ";
          }
          std::cout<<" }"<<std::endl;
        }
      }
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetPrintRecoClusterTruth);
