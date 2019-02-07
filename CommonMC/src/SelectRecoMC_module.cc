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
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "DataProducts/inc/IndexMap.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
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
      fhicl::Atom<art::InputTag> VDSPC { Name("VDSPCollection"),
	Comment("Virtual Detector StepPointMC collection")};
   };
    using Parameters = art::EDProducer::Table<Config>;
    explicit SelectRecoMC(const Parameters& conf);
    void produce(art::Event& evt) override;

  private:
  // utility functions
    void fillKalSeedMC(KalSeed const& seed, StrawDigiMCCollection const& sdmcc,
	StepPointMCCollection const& vdspc, PrimaryParticle const& pp, KalSeedMC& mcseed);
    int _debug;
    art::InputTag _pp, _ccc, _crvccc, _sdmcc, _crvdmcc, _vdspc;
    std::vector<std::string> _kff;
  };

  SelectRecoMC::SelectRecoMC(const Parameters& config )  : 
    _debug(config().debug()),
    _pp(config().PP()),
    _ccc(config().CCC()),
    _crvccc(config().CrvCCC()),
    _sdmcc(config().SDMCC()),
    _crvdmcc(config().CRVDMCC()),
    _vdspc(config().VDSPC()),
    _kff(config().KFFInstances())
  {
    consumes<PrimaryParticle>(_pp);
    consumesMany<KalSeedCollection>();
    consumes<CaloClusterCollection>(_ccc);
    consumes<CrvCoincidenceClusterCollection>(_crvccc);
    produces <IndexMap>("StrawDigiMCMap"); 
    produces <IndexMap>("CrvDigiMCMap"); 
    produces <KalSeedMCCollection>(); 
    produces <KalSeedMCAssns>();
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
    auto vdspch = event.getValidHandle<StepPointMCCollection>(_vdspc);
    auto const& vdspc = *vdspch;
    // some things needed for creating Ptrs before the collection is in the event
    auto KalSeedMCCollectionPID = getProductID<KalSeedMCCollection>();
    auto KalSeedMCCollectionGetter = event.productGetter(KalSeedMCCollectionPID);
// create output; these are pruned collections containing only
// products related to the reconstruction output or the event primary
    std::unique_ptr<IndexMap> sdmcim(new IndexMap);
    std::unique_ptr<IndexMap> crvdmcim(new IndexMap);
    std::unique_ptr<KalSeedMCCollection> ksmcc(new KalSeedMCCollection);
    std::unique_ptr<KalSeedMCAssns> ksmca(new KalSeedMCAssns);
    // straw hit (indices) that are referenced by the tracks
    std::set<StrawHitIndex> shindices;
    // CaloClusters that are referenced by the tracks
//    std::set<art::Ptr<CaloCluster> > caloclusters;
    // loop over input KalFinalFit products
    for (auto const& kff : _kff) {
    // get all products from this
      art::ModuleLabelSelector kffsel(kff);
      std::vector< art::Handle<KalSeedCollection> > seedhs;
      event.getMany<KalSeedCollection>(kffsel, seedhs);
      if(_debug > 1) std::cout << "Found " << seedhs.size() << " collections from module " << kff << std::endl;
      // loop over the KalSeeds and the hits inside them
      for(auto const& seedh : seedhs) {
	auto const& seedc = *seedh;
	if(_debug > 1) std::cout << "Found " << seedc.size() << " seeds from collection " << kff << std::endl;
	for(auto iseed=seedc.begin(); iseed!=seedc.end(); ++iseed){
	  auto const& seed = *iseed;
//	  if(seed.hasCaloCluster())caloclusters.insert(seed.caloCluster());
	  for( auto const& tsh : seed.hits() ) {
	    shindices.insert(tsh.index());
	  }
// create the KalSeedMC for this reco seed
	  KalSeedMC mcseed;
	  fillKalSeedMC(seed,sdmcc,vdspc,pp,mcseed);
	  ksmcc->push_back(mcseed);
	  // fill the Assns; this needs Ptrs
	  auto mcseedp = art::Ptr<KalSeedMC>(KalSeedMCCollectionPID,ksmcc->size()-1,KalSeedMCCollectionGetter);
	  auto seedp = art::Ptr<KalSeed>(seedh,std::distance(seedc.begin(),iseed));
	  ksmca->addSingle(seedp,mcseedp);
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
    event.put(std::move(ksmcc));
    event.put(std::move(ksmca));
  }

  void SelectRecoMC::fillKalSeedMC(KalSeed const& seed, 
    StrawDigiMCCollection const& sdmcc,
    StepPointMCCollection const& vdspc,
    PrimaryParticle const& pp,
    KalSeedMC& mcseed) {
  // find the associated SimParticles for this KalSeed
    std::vector<TrkMCTools::spcount> spcc;
    TrkMCTools::findMCTrk(seed,spcc,sdmcc);
    // find all the StepPointMCs associated with the primary particle
    if(spcc.size()>0){
    // create a stub for the primary particle
      auto const& pspc = spcc.front();
      SimPartStub pstub(pspc._spp);
      // these should be in the constructor: requires moving spcount declaration FIXME!
      pstub._nhits = pspc._count;
      pstub._nactive = pspc._acount;
      // no matching to true CaloCluster yet FIXME!
      // find the relationship with the primary particle (s).  Choose the closest
      for(auto const& spp : pp.primarySimParticles()) {
	MCRelationship prel(pspc._spp,spp);
	if(prel > pstub._rel) pstub._rel = prel;
      }
      mcseed._simps.push_back(pstub);
    // create stubs for the rest of the associated particles.  These refrence the primary
    // particle of this seed
      for(size_t isp=1;isp < spcc.size(); ++isp){
	auto const& spc = spcc[isp];
	SimPartStub stub(spc._spp);
	stub._nhits = spc._count;
	stub._nactive = spc._acount;
	stub._rel = MCRelationship(spc._spp,pspc._spp);
	mcseed._simps.push_back(pstub);
      }
// now find matching VD hits and record their info too
      for(auto const& vdsp : vdspc ) {
	if(vdsp.simParticle() == pspc._spp){
	  if(_debug > 1) std::cout << "Found matching VD StepPoint position" 
	  << vdsp.position() << " VDID = " << vdsp.virtualDetectorId() << std::endl;
	  mcseed._vdsteps.push_back(VDStep(vdsp));
	}
      }
    // find the SimParticle referenced by individual hits
      for(auto const& hit : seed.hits() ) {
	int spref(-1);
	auto const& sdmc = sdmcc.at(hit.index()); // bounds-check for security;
	for(size_t isp=1;isp < spcc.size(); ++isp){
	  auto const& spc = spcc[isp];
	  if(sdmc.earlyStepPointMC()->simParticle() == spc._spp){
	    spref = isp;
	    break;
	  }
	}
	// record the reference
	mcseed._digisimps.push_back(spref);
      }
    }
  }
}
DEFINE_ART_MODULE(mu2e::SelectRecoMC)
