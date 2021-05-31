//
//  Create a subset of reco-related MC objects (StrawDigiMC, etc) from
//  reconstruction output objects, as well as MC objects related
//  to the MC primary object.
//
// Original author: Dave Brown (LBNL) Feb 2019
// art
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "DataProducts/inc/IndexMap.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/CrvDigi.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "RecoDataProducts/inc/RecoCount.hh"
// Utilities
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

// C++
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
      fhicl::Atom<int>  debug                   { Name("debugLevel"),                     Comment("Debug Level"), 0};
      fhicl::Atom<bool> saveEnergySteps         { Name("SaveEnergySteps"),	          Comment("Save all StepPoints that contributed energy to a StrawDigi"), false};
      fhicl::Atom<bool> saveUnused              { Name("SaveUnusedDigiMCs"),	          Comment("Save StrawDigiMCs from particles used in any fit"), true};
      fhicl::Atom<bool> saveAllUnused           { Name("SaveAllUnusedDigiMCs"),	          Comment("Save all StrawDigiMCs from particles used in any fit"), false};
      fhicl::Atom<art::InputTag> PP             { Name("PrimaryParticle"),	          Comment("PrimaryParticle")};
      fhicl::Atom<art::InputTag> CCC            { Name("CaloClusterCollection"),          Comment("CaloClusterCollection")};
      fhicl::Atom<art::InputTag> CrvCCC         { Name("CrvCoincidenceClusterCollection"),Comment("CrvCoincidenceClusterCollection")};
      fhicl::Atom<art::InputTag> SDC            { Name("StrawDigiCollection"),	          Comment("StrawDigiCollection")};
      fhicl::Atom<art::InputTag> SHFC           { Name("StrawHitFlagCollection"),         Comment("StrawHitFlagCollection")};
      fhicl::Atom<art::InputTag> CHC            { Name("ComboHitCollection"),	          Comment("ComboHitCollection for the original StrawHits (not Panel hits)")};
      fhicl::Atom<art::InputTag> CDC            { Name("CaloDigiCollection"),	          Comment("CaloDigiCollection")};
      fhicl::Atom<art::InputTag> SDMCC          { Name("StrawDigiMCCollection"),          Comment("StrawDigiMCCollection")};
      fhicl::Atom<art::InputTag> CRVDC          { Name("CrvDigiCollection"),	          Comment("CrvDigiCollection")};
      fhicl::Atom<art::InputTag> CRVDMCC        { Name("CrvDigiMCCollection"),	          Comment("CrvDigiMCCollection")};
      fhicl::Sequence<std::string> KalSeeds     { Name("KalSeedCollections"),	          Comment("KalSeedCollections")};
      fhicl::Sequence<std::string> HelixSeeds   { Name("HelixSeedCollections"),	          Comment("HelixSeedCollections")};
      fhicl::Atom<art::InputTag> VDSPC          { Name("VDSPCollection"),	          Comment("Virtual Detector StepPointMC collection")};
      fhicl::Sequence<art::InputTag> SPTO       { Name("TimeOffsets"),	                  Comment("Sim Particle Time Offset Maps")};
      fhicl::Atom<double> CCME                  { Name("CaloClusterMinE"),                Comment("Minimum energy CaloCluster to save digis (MeV)")};
   };
   using Parameters = art::EDProducer::Table<Config>;
   explicit SelectRecoMC(const Parameters& conf);
   void produce(art::Event& evt) override;
        
  private:
   typedef std::vector<TrkMCTools::spcount>SPCC;
   typedef std::set<StrawHitIndex> SHIS;
  // utility functions
   void fillTSHMC         (KalSeed const& seed, SPCC const& spcc, StrawDigiMCCollection const& sdmcc,KalSeedMC& mcseed);
   void fillUnusedTSHMC   (SPCC const& spcc, StrawDigiMCCollection const& sdmcc, KalSeedMC& mcseed);
   void fillVDSP          (GeomHandle<DetectorSystem>const& det, art::Ptr<SimParticle> const& psp, StepPointMCCollection const& vdspc, KalSeedMC& mcseed);
   void fillSPStubs       (SPCC const& spcc, PrimaryParticle const& pp, KalSeedMC& mcseed);
   void fillSDMCI         (KalSeedMC const& mcseed,SHIS& shindices);
   void fillStrawHitCounts(ComboHitCollection const& chc, StrawHitFlagCollection const& shfc, RecoCount& nrec);
   void fillTrk           (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, PrimaryParticle const& pp, RecoCount& nrec);
   void fillCrv           (art::Event& event, PrimaryParticle const& pp, RecoCount& nrec);
   void fillCalo          (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs,PrimaryParticle const& pp, RecoCount& nrec);
   int _debug;
   bool _saveallenergy, _saveunused, _saveallunused;
   art::InputTag _pp, _ccc, _crvccc, _sdc, _shfc, _chc, _cdc, _sdmcc, _crvdc, _crvdmcc, _vdspc;
   std::vector<std::string> _kscs, _hscs;
   SimParticleTimeOffset _toff;
   double _ccme;
   // cache
   double _mbtime; // period of 1 microbunch

  };

  SelectRecoMC::SelectRecoMC(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _saveallenergy(config().saveEnergySteps()),
    _saveunused(config().saveUnused()),
    _saveallunused(config().saveAllUnused()),
    _pp(config().PP()),
    _ccc(config().CCC()),
    _crvccc(config().CrvCCC()),
    _sdc(config().SDC()),
    _shfc(config().SHFC()),
    _chc(config().CHC()),
    _cdc(config().CDC()),
    _sdmcc(config().SDMCC()),
    _crvdc(config().CRVDC()),
    _crvdmcc(config().CRVDMCC()),
    _vdspc(config().VDSPC()),
    _kscs(config().KalSeeds()),
    _hscs(config().HelixSeeds()),
    _toff(config().SPTO()),
    _ccme(config().CCME())
  {
    consumes<StrawDigiCollection>(_sdc);
    consumes<StrawDigiADCWaveformCollection>(_sdc);
    consumes<StrawHitFlagCollection>(_shfc);
    consumes<ComboHitCollection>(_chc);
    consumes<CaloDigiCollection>(_cdc);
    consumes<CaloClusterCollection>(_ccc);
    consumes<CrvDigiCollection>(_crvdc);
    consumesMany<KalSeedCollection>();
    consumes<CrvCoincidenceClusterCollection>(_crvccc);
    consumes<PrimaryParticle>(_pp);
    consumes<StrawDigiMCCollection>(_sdmcc);
    consumes<CrvDigiMCCollection>(_crvdmcc);
    produces <IndexMap>("StrawDigiMap"); 
    produces <IndexMap>("CrvDigiMap"); 
    produces <KalSeedMCCollection>(); 
    produces <KalSeedMCAssns>();    
    produces <CaloDigiCollection>();    
    produces <StrawDigiCollection>();
    produces <StrawDigiADCWaveformCollection>();
    produces <StrawHitFlagCollection>();
    produces <CrvDigiCollection>();
    produces <CrvRecoPulseCollection>();
    produces <CrvCoincidenceClusterCollection>();
    produces <RecoCount>();
    
    if (_debug > 0)
    {
       std::cout << "Using KalSeed collections from ";
       for (auto const& kff : _kscs) std::cout << kff << " " << std::endl;
      
       std::cout << "Using HelixSeed collections from ";
       for (auto const& hsc : _hscs) std::cout << hsc << " " << std::endl;
    }
  }

  void SelectRecoMC::produce(art::Event& event) {
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);

    auto pph = event.getValidHandle<PrimaryParticle>(_pp);
    auto const& pp = *pph;

    std::unique_ptr<RecoCount> nrec(new RecoCount);
    std::set<art::Ptr<CaloCluster> > ccptrs;
    fillTrk(event,ccptrs,pp,*nrec.get());
    fillCrv(event, pp, *nrec.get());
    fillCalo(event, ccptrs, pp, *nrec.get());

    event.put(std::move(nrec));
  }

  void SelectRecoMC::fillSPStubs(SPCC const& spcc, PrimaryParticle const& pp, KalSeedMC& mcseed) {
    // find all the StepPointMCs associated with the primary particle
    if(spcc.size()>0){
    // create a stub for the primary particle
      auto const& pspc = spcc.front();
      SimPartStub pstub(pspc._spp);
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
    // particle of this seed for the relationship
      for(size_t isp=1;isp < spcc.size(); ++isp){
	auto const& spc = spcc[isp];
	SimPartStub stub(spc._spp);
	stub._nhits = spc._count;
	stub._nactive = spc._acount;
	stub._rel = MCRelationship(spc._spp,pspc._spp);
	mcseed._simps.push_back(stub);
      }
    }
  }
  
  void SelectRecoMC::fillTSHMC(KalSeed const& seed, SPCC const& spcc,
      StrawDigiMCCollection const& sdmcc, KalSeedMC& mcseed) {
    for(auto const& hit : seed.hits() ) {
      // create a TrkStrawHitMC for each hit on the seed
      TrkStrawHitMC tshmc;
      tshmc._sdmcindex = hit.index();
      // find the referenced sim particle
      int spref(-1);
      auto const& sdmc = sdmcc.at(hit.index()); // bounds-check for security;
      for(size_t isp=0;isp < spcc.size(); isp++){
	auto const& spc = spcc[isp];
	if(sdmc.earlyStrawGasStep()->simParticle() == spc._spp){
	  spref = isp;
	  break;
	}
      }
      if(spref < 0)throw cet::exception("Reco")<<"mu2e::SelectRecoMC: missing index"<< std::endl;
      tshmc._spindex = spref;
      // fill other info directly from the StrawDigiMC
      tshmc._energySum = sdmc.triggerEnergySum(sdmc.earlyEnd());
      const auto& mcstep = *(sdmc.earlyStrawGasStep());
      tshmc._cpos = Geom::toXYZVec(sdmc.clusterPosition(sdmc.earlyEnd()));
      tshmc._mom = mcstep.momentum();
      tshmc._time = fmod(mcstep.time(),_mbtime);
      tshmc._strawId = sdmc.strawId();
      mcseed._tshmcs.push_back(tshmc);
    }
  }

  void SelectRecoMC::fillUnusedTSHMC( SPCC const& spcc,
      StrawDigiMCCollection const& sdmcc,
      KalSeedMC& mcseed) {
    // either keep hits only from the primary or from all contributing
    size_t ispmax = _saveallunused? spcc.size() : 1;
    for(size_t isp=0; isp < ispmax; ++isp){
      auto const& spc = spcc[isp];
      for (size_t isdmc=0; isdmc < sdmcc.size(); isdmc++){
	auto const& sdmc = sdmcc[isdmc];
	auto const& mcstep = *(sdmc.earlyStrawGasStep());
	if(mcstep.simParticle() == spc._spp){
	  // search to see if the associated digi is already on the track
	  bool used(false);
	  for(auto const& tshmc : mcseed._tshmcs ) {
	    if(isdmc == tshmc.strawDigiMCIndex()){
	      used = true;
	      break;
	    }
	  }
	  if(!used){
	    // record the reference
	    TrkStrawHitMC tshmc;
	    tshmc._sdmcindex = isdmc;
	    tshmc._spindex = isp;
	    tshmc._energySum = sdmc.triggerEnergySum(sdmc.earlyEnd());
	    tshmc._cpos = Geom::toXYZVec(sdmc.clusterPosition(sdmc.earlyEnd()));
	    tshmc._mom = mcstep.momentum();
	    tshmc._time = sdmc.clusterTime(sdmc.earlyEnd());
	    tshmc._strawId = sdmc.strawId();
	    mcseed._tshmcs.push_back(tshmc);
	  }
	}
      }
    }
  }

  void SelectRecoMC::fillSDMCI(KalSeedMC const& mcseed, SHIS& shindices) {
    for(auto tshmc : mcseed.trkStrawHitMCs()) {
      shindices.insert(tshmc.strawDigiMCIndex());
    }
  }

  void SelectRecoMC::fillVDSP( GeomHandle<DetectorSystem>const& det,
      art::Ptr<SimParticle> const& psp,
      StepPointMCCollection const& vdspc, KalSeedMC& mcseed) {
    // loop over all the StepPointMCs and pick the ones that have his
    // for the primary particle
    for(auto const& vdsp : vdspc ) {
      if(vdsp.simParticle() == psp){
	if(_debug > 1) std::cout << "Found matching VD StepPoint position" 
	  << vdsp.position() << " momentum " << vdsp.momentum() 
	    << " time " << vdsp.time() << " VDID = " << vdsp.virtualDetectorId() << std::endl;
	VDStep vds(det->toDetector(vdsp.position()),
	    vdsp.momentum() ,
	    _toff.timeWithOffsetsApplied(vdsp),
	    vdsp.virtualDetectorId());
	mcseed._vdsteps.push_back(vds);
      }
    }
  }

  void SelectRecoMC::fillStrawHitCounts(ComboHitCollection const& chc, StrawHitFlagCollection const& shfc, RecoCount& nrec) {
// test
    if(chc.size() != shfc.size())
      throw cet::exception("Reco")<<"mu2e::SelectRecoMC: inconsistent collections"<< std::endl; 
    for(const auto& shf : shfc) {
      if(shf.hasAllProperties(StrawHitFlag::energysel))++nrec._nshfesel;
      if(shf.hasAllProperties(StrawHitFlag::radsel))++nrec._nshfrsel;
      if(shf.hasAllProperties(StrawHitFlag::timesel))++nrec._nshftsel;
      if(shf.hasAllProperties(StrawHitFlag::bkg))++nrec._nshfbkg;
      if(shf.hasAllProperties(StrawHitFlag::trksel))++nrec._nshftpk;
    }
    nrec._nstrawdigi = chc.size();
    // fill straw hit time histogram
    for(auto const& ch : chc)nrec._shthist.fill(ch.time());
  }

  void SelectRecoMC::fillTrk( art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs,
      PrimaryParticle const& pp, RecoCount& nrec) {
    GeomHandle<DetectorSystem> det;
     // Tracker-reated data products
    auto sdch = event.getValidHandle<StrawDigiCollection>(_sdc);
    auto const& sdc = *sdch;
    auto sdadcch = event.getValidHandle<StrawDigiADCWaveformCollection>(_sdc);
    auto const& sdadcc = *sdadcch;
    auto shfch = event.getValidHandle<StrawHitFlagCollection>(_shfc);
    auto const& shfc = *shfch;
    auto chch = event.getValidHandle<ComboHitCollection>(_chc);
    auto const& chc = *chch;
    auto sdmcch = event.getValidHandle<StrawDigiMCCollection>(_sdmcc);
    auto const& sdmcc = *sdmcch;
    // virtual detector hits won't exist in no-primary samples: don't flag this as an error
    art::Handle<StepPointMCCollection> vdspch;
    event.getByLabel<StepPointMCCollection>(_vdspc,vdspch);
    // some things needed for creating Ptrs before the collection is in the event
    auto KalSeedMCCollectionPID = event.getProductID<KalSeedMCCollection>();
    auto KalSeedMCCollectionGetter = event.productGetter(KalSeedMCCollectionPID);
// create products related to the reconstruction output or the event primary
    std::unique_ptr<KalSeedMCCollection> ksmcc(new KalSeedMCCollection);
    std::unique_ptr<KalSeedMCAssns> ksmca(new KalSeedMCAssns);
    std::unique_ptr<StrawDigiCollection> ssdc(new StrawDigiCollection);
    std::unique_ptr<StrawDigiADCWaveformCollection> ssdadcc(new StrawDigiADCWaveformCollection);
    std::unique_ptr<StrawHitFlagCollection> sshfc(new StrawHitFlagCollection);
// index maps between original collections and pruned collections
    std::unique_ptr<IndexMap> sdmcim(new IndexMap);
    // straw hit (indices) that are referenced by the tracks, or particles
    // that contributed to the track
    SHIS shindices;
    // set of SimParticles to save for this event
    std::set<art::Ptr<SimParticle> > simps;
    // add the MC primary SimParticles
    for(auto const& spp : pp.primarySimParticles()) simps.insert(spp);
    // loop over input KalFinalFit products
    for (auto const& ksc : _kscs) {
    // get all products from this
      art::ModuleLabelSelector kscsel(ksc);
      std::vector< art::Handle<KalSeedCollection> > seedhs = event.getMany<KalSeedCollection>(kscsel);
      if(_debug > 1) std::cout << "Found " << seedhs.size() << " collections from module " << ksc << std::endl;
      // loop over the KalSeeds and the hits inside them
      for(auto const& seedh : seedhs) {
	auto const& seedc = *seedh;
	if(_debug > 1) std::cout << "Found " << seedc.size() << " seeds from collection " << ksc << std::endl;
	for(auto iseed=seedc.begin(); iseed!=seedc.end(); ++iseed){
	  auto const& seed = *iseed;
	  // find the associated SimParticles for this KalSeed.  This ranks them by their hit count
	  SPCC spcc;
	  TrkMCTools::findMCTrk(seed,spcc,sdmcc,_saveallenergy);
	  if(_debug > 2) std::cout << "Found " << spcc.size() << " Associated SimParticles for KalSeed " << std::endl;
	  // add these to the set of SimParticles to save
	  for(auto const& spc : spcc) simps.insert(spc._spp);
	  // create the KalSeedMC for this reco seed and fill the parts
	  KalSeedMC mcseed;
	  fillSPStubs(spcc,pp,mcseed);
	  fillTSHMC(seed,spcc,sdmcc,mcseed);
	  // add DigiMCs not used in the track but from true particles used in the track
	  if(_saveunused)fillUnusedTSHMC(spcc,sdmcc,mcseed);
	  if(spcc.size() > 0 && vdspch.isValid()){
	    auto const& vdspc = *vdspch;
	    fillVDSP(det,spcc.front()._spp,vdspc,mcseed);
	  }
	  ksmcc->push_back(mcseed);
	  // fill indices from all digis; those on the track and those from the MC true particle too
	  fillSDMCI(mcseed,shindices);
	  // fill the Assns; this needs Ptrs
	  auto mcseedp = art::Ptr<KalSeedMC>(KalSeedMCCollectionPID,ksmcc->size()-1,KalSeedMCCollectionGetter);
	  auto seedp = art::Ptr<KalSeed>(seedh,std::distance(seedc.begin(),iseed));
	  ksmca->addSingle(seedp,mcseedp);
	  // record the CaloCluster associated with this seed (if any)
	  if(seed.hasCaloCluster())ccptrs.insert(seed.caloCluster());
	  if(_debug > 2) std::cout << "KalSeedMC has " << mcseed._tshmcs.size() 
	    << " MCHits , KalSeed has " << seed.hits().size() << " TrkStrawHits" << std::endl;
	} 
      }
    }
    // get straw indices from all helices too
    for (auto const& hsc : _hscs) {
    // get all products from this
      art::ModuleLabelSelector hscsel(hsc);
      std::vector< art::Handle<HelixSeedCollection> > seedhs = event.getMany<HelixSeedCollection>(hscsel);
      if(_debug > 1) std::cout << "Found " << seedhs.size() << " collections from module " << hsc << std::endl;
      // loop over the HelixSeeds and the hits inside them
      for(auto const& seedh : seedhs) {
	auto const& seedc = *seedh;
	if(_debug > 1) std::cout << "Found " << seedc.size() << " seeds from collection " << hsc << std::endl;
	for(auto iseed=seedc.begin(); iseed!=seedc.end(); ++iseed){
	  auto const& seed = *iseed;
	  // go back to StrawHit indices (== digi indices for reco)
	  std::vector<StrawHitIndex> shids;
	  for(size_t ihit = 0; ihit < seed.hits().size(); ihit++)
	    seed.hits().fillStrawHitIndices(event,ihit,shids);
	  // add these to the set (duplicates are suppressed)
	  for(auto shid : shids)
	    shindices.insert(shid);
	}
      }
    }
  // fill the StrawIndex map with the complete list of indices.
    StrawHitIndex shcount(0);
    ssdc->reserve(shindices.size());
    ssdadcc->reserve(shindices.size());
    for(auto shindex : shindices){
      sdmcim->addElement(shindex,shcount++);
// deep-copy the selected StrawDigis and StrawHitFlags
      ssdc->push_back(sdc[shindex]);
      ssdadcc->push_back(sdadcc[shindex]);
      sshfc->push_back(shfc[shindex]);
    }
    if(_debug > 1) std::cout << "Selected " << shcount << " StrawDigis" << std::endl;

    // fill detailed StrawHit counts
    fillStrawHitCounts(chc,shfc,nrec);
    event.put(std::move(sdmcim),"StrawDigiMap");
    event.put(std::move(ssdc));
    event.put(std::move(ssdadcc));
    event.put(std::move(sshfc));
    event.put(std::move(ksmcc));
    event.put(std::move(ksmca));
  }

  void SelectRecoMC::fillCrv(art::Event& event, 
      PrimaryParticle const& pp, RecoCount& nrec) {
// find Crv data in event
    auto crvccch = event.getValidHandle<CrvCoincidenceClusterCollection>(_crvccc);
    auto const& crvccc = *crvccch;
    auto crvdch = event.getValidHandle<CrvDigiCollection>(_crvdc);
    auto const& crvdc = *crvdch;
    auto crvdmcch = event.getValidHandle<CrvDigiMCCollection>(_crvdmcc);
    auto const& crvdmcc = *crvdmcch;
    auto CrvRecoPulseCollectionPID = event.getProductID<CrvRecoPulseCollection>();
    auto CrvRecoPulseCollectionGetter = event.productGetter(CrvRecoPulseCollectionPID);
// create new Crv collections
    std::unique_ptr<CrvDigiCollection> scrvdc(new CrvDigiCollection);
    std::unique_ptr<CrvRecoPulseCollection> scrvrpc(new CrvRecoPulseCollection);
    std::unique_ptr<CrvCoincidenceClusterCollection> scrvccc(new CrvCoincidenceClusterCollection);
    std::unique_ptr<IndexMap> crvdmcim(new IndexMap);
// loop over CrvCoincidenceClusters
    std::set<uint16_t> crvindices;
    for(auto const& crvcc: crvccc) {
      std::vector<art::Ptr<CrvRecoPulse>> pulses;
      for(auto const& crvrp : crvcc.GetCrvRecoPulses()){
      // deep-copy the pulses used in coincidences: we do NOT update the digi indies,
      // the map must be used to connect them
	scrvrpc->push_back(*crvrp);
	auto crvrpp = art::Ptr<CrvRecoPulse>(CrvRecoPulseCollectionPID,scrvrpc->size()-1,CrvRecoPulseCollectionGetter);
	pulses.push_back(crvrpp);
	for(auto index : crvrp->GetWaveformIndices()){
	  crvindices.insert(index);
	}
      }
    // deep-copy the coincidence-cluster with updated Reco Pulses
      CrvCoincidenceCluster scrvcc(crvcc);
      scrvcc.SetCrvRecoPulses(pulses);
      scrvccc->push_back(scrvcc);
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
    for(auto crvindex : crvindices){
      crvdmcim->addElement(crvindex,crvcount++);
      // deep-copy the selected CrvDigis
      scrvdc->push_back(crvdc.at(crvindex));
    }
    // update reco count
    nrec._ncrvdigi = crvdc.size();
    // put new data in event
    event.put(std::move(scrvdc));
    event.put(std::move(scrvrpc));
    event.put(std::move(scrvccc));
    event.put(std::move(crvdmcim),"CrvDigiMap");
  }

  void SelectRecoMC::fillCalo(art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, 
                              PrimaryParticle const& pp, RecoCount& nrec)
  {

    auto cdch = event.getValidHandle<CaloDigiCollection>(_cdc);
    auto const& cdc = *cdch;
    auto ccch = event.getValidHandle<CaloClusterCollection>(_ccc);
    auto const& ccc = *ccch;
    
    std::unique_ptr<CaloDigiCollection> scdc(new CaloDigiCollection);
    // reco count
    nrec._ncalodigi = cdc.size();
    nrec._ncc = ccc.size();
    nrec._cce = 0.0;
    // loop over all the CaloClusters and mark the ones that are above energy for saving by adding their Ptrs to the list
    for(unsigned icc=0;icc < ccc.size(); icc++){
      auto const& cc = ccc[icc];
      nrec._cce += cc.energyDep();
      if(cc.energyDep() > _ccme){
	auto ccp = art::Ptr<CaloCluster>(ccch,icc);
	ccptrs.insert(ccp);
      }
    }
    // deep-copy CaloDigis from selected clusters
    for(auto const& ccptr : ccptrs) {
      for(auto const& cchptr : ccptr->caloHitsPtrVector()){
	for (auto const& rcdptr : cchptr->recoCaloDigis()){
          // deep-copy CaloDigis used in clusters
	  scdc->push_back(*rcdptr->caloDigiPtr());
	}
      }
    }
    // put new data into event
    event.put(std::move(scdc));
  }

}
DEFINE_ART_MODULE(mu2e::SelectRecoMC)
