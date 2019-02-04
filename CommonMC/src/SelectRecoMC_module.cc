//
//  Create a subset of reco-related MC objects (StrawDigiMC, etc) from
//  reconstruction output objects, as well as MC objects related
//  to the MC primary object.
//
// Original author: David Brown (LBNL) Feb 2019
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e 
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "DataProducts/inc/IndexMap.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <string>

namespace mu2e {
  class SelectRecoMC : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> debug{ Name("debugLevel"),
	Comment("Debug Level"), 0};
      fhicl::Atom<art::InputTag> PP { Name("PrimaryParticle"),
	Comment("PrimaryParticle producer")};
      fhicl::Atom<art::InputTag> CCC { Name("CaloClusterCollection"),
	Comment("CaloClusterCollection producer")};
      fhicl::Atom<art::InputTag> CrvCCC { Name("CrvCoincidenceClusterCollection"),
	Comment("CrvCoincidenceClusterCollection producer")};
      fhicl::Atom<art::InputTag> SDMCC { Name("StrawDigiMCCollection"),
	Comment("StrawDigiMCCollection producer")};
      fhicl::Atom<art::InputTag> CRVDMCC { Name("CrvDigiMCCollection"),
	Comment("CrvDigiMCCollection producer")};
      fhicl::Sequence<std::string> KFFInstances { Name("KFFInstances"),
	Comment("KalFinalFit Module Instances")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit SelectRecoMC(const Parameters& conf);
    void produce(art::Event& evt) override;

  private:
    int _debug;
    art::InputTag _pp, _ccc, _crvccc, _sdmcc, _crvdmcc;
    std::vector<std::string> _kff; 
  };

  SelectRecoMC::SelectRecoMC(const Parameters& config )  : 
    _debug(config().debug()),
    _pp(config().PP()),
    _ccc(config().CCC()),
    _crvccc(config().CrvCCC()),
    _sdmcc(config().SDMCC()),
    _crvdmcc(config().CRVDMCC())
  {
    consumes<PrimaryParticle>(_pp);
    consumesMany<KalSeedCollection>();
    consumes<CaloClusterCollection>(_ccc);
    consumes<CrvCoincidenceClusterCollection>(_crvccc);
    produces <IndexMap>("StrawDigiMCMap"); 
    produces <IndexMap>("CrvDigiMCMap"); 
    if(_debug > 0){
      std::cout << "Using KalSeed collections from ";
      for (auto const& kff : _kff)
	std::cout << kff << " " << std::endl;
    }
  }


  void SelectRecoMC::produce(art::Event& event) {
//  Find the inputs: the MC primary object
    auto pph = event.getValidHandle<PrimaryParticle>(_pp);
    auto const& pp = *pph;
    auto sdmcch = event.getValidHandle<StrawDigiMCCollection>(_sdmcc);
    auto const& sdmcc = *sdmcch;
    auto crvccch = event.getValidHandle<CrvCoincidenceClusterCollection>(_crvccc);
    auto const& crvccc = *crvccch;
    auto crvdmcch = event.getValidHandle<CrvDigiMCCollection>(_crvdmcc);
    auto const& crvdmcc = *crvdmcch;
// create output; these are pruned collections containing only
// products related to the reconstruction output or the event primary
    std::unique_ptr<IndexMap> sdmcim(new IndexMap);
    std::unique_ptr<IndexMap> crvdmcim(new IndexMap);
    // straw hit (indices) that are referenced by the tracks
    std::set<StrawHitIndex> shindices;
    // CaloClusters that are referenced by the tracks
    std::set<art::Ptr<CaloCluster> > caloclusters;
    // loop over input KalFinalFit products
    for (auto const& kff : _kff) {
    // get all products from this
      art::ModuleLabelSelector kffsel(kff);
      std::vector< art::Handle<KalSeedCollection> > seedhs;
      event.getMany<KalSeedCollection>(kffsel, seedhs);
      // loop over the KalSeeds and the hits inside them
      for(auto const& seedh : seedhs) {
	auto const& seedc = *seedh;
	for(auto const& seed : seedc) {
	  if(seed.hasCaloCluster())caloclusters.insert(seed.caloCluster());
	  for( auto const& tsh : seed.hits() ) {
	    shindices.insert(tsh.index());
	  }
	} 
      }
    }
    // From these, create the vector of SimParticles
    std::set<art::Ptr<SimParticle> > simps;
    for(auto shindex : shindices) {
      auto const& sdmc = sdmcc[shindex];
      // look at all contributing StepPoints.  Might be overkill, but
      // guarantees all possible functionality
      for(auto const& spmc : sdmc.stepPointMCs()){
	simps.insert(spmc->simParticle());
      }
    }
    // add the MC primary SimParticles
    for(auto const& spp : pp.primarySimParticles())
      simps.insert(spp);
    // Should add the SimParticles from the CaloClusters and CRV Coincidences FIXME!!
    //
    //  Find all the StrawDigiMCs that reference these SimParticles and add them to the index list
    //  This guarantees all true digis from the relevant particles are saved, even if they weren't
    //  used in the track fit
    for(auto sdmci=sdmcc.begin();sdmci!=sdmcc.end(); ++sdmci) {
      auto const& sdmc = *sdmci;
      for(auto const& spmc : sdmc.stepPointMCs()){
	if(std::find(simps.begin(),simps.end(),spmc->simParticle()) != simps.end()){
	  shindices.insert(std::distance(sdmcc.begin(),sdmci));
	}
      }
    }
    // fill the StrawIndex map with the complete list of indices.  Note that std::set has kept these
    // sorted for us, nice!
    StrawHitIndex shcount(0);
    for(auto shindex : shindices) sdmcim->addElement(shindex,shcount++);
    // now the same thing for the CRV
    std::set<uint16_t> crvindices;
    for(auto const& crvcc: crvccc) {
      for(auto const& crvrp : crvcc.GetCrvRecoPulses()){
	for(auto index : crvrp->GetWaveformIndices()){
	  crvindices.insert(index);
	}
      }
    }
    // add indices for digis associated with the MC primary particle(s)
    for(auto icrv = crvdmcc.begin(); icrv != crvdmcc.end();++icrv) {
      auto const& crvdmc = *icrv;
      if(std::find(pp.primarySimParticles().begin(),
	    pp.primarySimParticles().end(),
	    crvdmc.GetSimParticle()) !=  pp.primarySimParticles().end()){
	crvindices.insert(std::distance(crvdmcc.begin(),icrv));
      }
    }
    // Fill CrvIndex map
    uint16_t crvcount(0);
    for(auto crvindex : crvindices) crvdmcim->addElement(crvindex,crvcount++);
    // do something with the CaloClusters: FIXME!

    // put output in event
    event.put(std::move(sdmcim),"StrawDigiMCMap");
    event.put(std::move(crvdmcim),"CrvDigiMCMap");
  }
}
DEFINE_ART_MODULE(mu2e::SelectRecoMC)
