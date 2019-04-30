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
#include "MCDataProducts/inc/StepPointMC.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/CrvDigi.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
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
      fhicl::Atom<int> debug{ Name("debugLevel"),
	Comment("Debug Level"), 0};
      fhicl::Atom<bool> saveEnergySteps{ Name("SaveEnergySteps"),
	Comment("Save all StepPoints that contributed energy to a StrawDigi"), false};
      fhicl::Atom<bool> saveUnused{ Name("SaveUnusedDigiMCs"),
	Comment("Save StrawDigiMCs from particles used in any fit"), true};
      fhicl::Atom<bool> saveAllUnused{ Name("SaveAllUnusedDigiMCs"),
	Comment("Save all StrawDigiMCs from particles used in any fit"), false};
      fhicl::Atom<art::InputTag> PP { Name("PrimaryParticle"),
	Comment("PrimaryParticle producer")};
      fhicl::Atom<art::InputTag> CCC { Name("CaloClusterCollection"),
	Comment("CaloClusterCollection producer")};
      fhicl::Atom<art::InputTag> CrvCCC { Name("CrvCoincidenceClusterCollection"),
	Comment("CrvCoincidenceClusterCollection producer")};
      fhicl::Atom<art::InputTag> SDC { Name("StrawDigiCollection"),
	Comment("StrawDigiCollection producer")};
      fhicl::Atom<art::InputTag> SHFC { Name("StrawHitFlagCollection"),
	Comment("StrawHitFlagCollection producer")};
      fhicl::Atom<art::InputTag> CDC { Name("CaloDigiCollection"),
	Comment("CaloDigiCollection producer")};
       fhicl::Atom<art::InputTag> SDMCC { Name("StrawDigiMCCollection"),
	Comment("StrawDigiMCCollection producer")};
      fhicl::Atom<art::InputTag> CRVDC { Name("CrvDigiCollection"),
	Comment("CrvDigiCollection producer")};
      fhicl::Atom<art::InputTag> CRVDMCC { Name("CrvDigiMCCollection"),
	Comment("CrvDigiMCCollection producer")};
      fhicl::Sequence<std::string> KFFInstances { Name("KFFInstances"),
	Comment("KalFinalFit Module Instances")};
      fhicl::Atom<art::InputTag> VDSPC { Name("VDSPCollection"),
	Comment("Virtual Detector StepPointMC collection")};
      fhicl::Sequence<art::InputTag> SPTO { Name("TimeOffsets"),
	Comment("Sim Particle Time Offset Maps")};
      fhicl::Atom<art::InputTag> CSSC { Name("CSSCollection"),
	Comment("CaloShowerSim collection")};
      fhicl::Atom<double> CCMCDT { Name("CaloClusterMCDTime"),
	Comment("Max time difference between CaloCluster and CaloShowerSim for MC match")};
      fhicl::Atom<double> CCMCE { Name("CaloMinE"),
	Comment("Min CaloShowerSim MC energy to include")};
 
   };
   using Parameters = art::EDProducer::Table<Config>;
   explicit SelectRecoMC(const Parameters& conf);
   void produce(art::Event& evt) override;

  private:
   typedef std::vector<TrkMCTools::spcount>SPCC;
   typedef std::set<StrawHitIndex> SHIS;
  // utility functions
   void fillTSHMC(KalSeed const& seed, SPCC const& spcc, StrawDigiMCCollection const& sdmcc,
        KalSeedMC& mcseed);
   void fillUnusedTSHMC( SPCC const& spcc, StrawDigiMCCollection const& sdmcc, KalSeedMC& mcseed);
   void fillVDSP( GeomHandle<DetectorSystem>const& det, art::Ptr<SimParticle> const& psp,
       StepPointMCCollection const& vdspc, KalSeedMC& mcseed);
   void fillSPStubs(SPCC const& spcc, PrimaryParticle const& pp, KalSeedMC& mcseed);
   void fillSDMCI(KalSeedMC const& mcseed,SHIS& shindices);
   void fillCaloClusterMC(CaloCluster const& cc, CaloShowerSimCollection const& cssc,
       PrimaryParticle const& pp, CaloClusterMC& ccmc);
   void fillStrawHitCounts(StrawHitFlagCollection const& shfc, RecoCount& nrec);
   int _debug;
   bool _saveallenergy, _saveunused, _saveallunused;
   art::InputTag _pp, _ccc, _crvccc, _sdc, _shfc, _cdc, _sdmcc, _crvdc, _crvdmcc, _vdspc, _cssc;
   std::vector<std::string> _kff;
   SimParticleTimeOffset _toff;
   double _ccmcdt, _ccmce;
   // cache
   double _mbtime; // period of 1 microbunch

  };

  SelectRecoMC::SelectRecoMC(const Parameters& config )  : 
    _debug(config().debug()),
    _saveallenergy(config().saveEnergySteps()),
    _saveunused(config().saveUnused()),
    _saveallunused(config().saveAllUnused()),
    _pp(config().PP()),
    _ccc(config().CCC()),
    _crvccc(config().CrvCCC()),
    _sdc(config().SDC()),
    _shfc(config().SHFC()),
    _cdc(config().CDC()),
    _sdmcc(config().SDMCC()),
    _crvdc(config().CRVDC()),
    _crvdmcc(config().CRVDMCC()),
    _vdspc(config().VDSPC()),
    _cssc(config().CSSC()),
    _kff(config().KFFInstances()),
    _toff(config().SPTO()),
    _ccmcdt(config().CCMCDT()),
    _ccmce(config().CCMCE())
  {
    consumes<StrawDigiCollection>(_sdc);
    consumes<StrawHitFlagCollection>(_shfc);
    consumes<CaloDigiCollection>(_cdc);
    consumes<CrvDigiCollection>(_crvdc);
    consumesMany<KalSeedCollection>();
    consumes<CaloClusterCollection>(_ccc);
    consumes<CrvCoincidenceClusterCollection>(_crvccc);
    consumes<PrimaryParticle>(_pp);
    consumes<StrawDigiMCCollection>(_sdmcc);
    consumes<CrvDigiMCCollection>(_crvdmcc);
    produces <IndexMap>("StrawDigiMap"); 
    produces <IndexMap>("CrvDigiMap"); 
    produces <IndexMap>("CaloDigiMap"); 
    produces <KalSeedMCCollection>(); 
    produces <KalSeedMCAssns>();
    produces <CaloClusterMCCollection>(); 
    produces <CaloClusterMCAssns>();
    produces <StrawDigiCollection>();
    produces <StrawHitFlagCollection>();
    produces <CaloDigiCollection>();
    produces <CrvDigiCollection>();
    produces <RecoCount>();
    if(_debug > 0){
      std::cout << "Using KalSeed collections from ";
      for (auto const& kff : _kff)
	std::cout << kff << " " << std::endl;
    }
  }

  void SelectRecoMC::produce(art::Event& event) {
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
    GeomHandle<DetectorSystem> det;
// update the time maps
    _toff.updateMap(event);
//  Find the inputs: the MC primary object
    auto pph = event.getValidHandle<PrimaryParticle>(_pp);
    auto const& pp = *pph;
    auto sdch = event.getValidHandle<StrawDigiCollection>(_sdc);
    auto const& sdc = *sdch;
    auto shfch = event.getValidHandle<StrawHitFlagCollection>(_shfc);
    auto const& shfc = *shfch;
    auto sdmcch = event.getValidHandle<StrawDigiMCCollection>(_sdmcc);
    auto const& sdmcc = *sdmcch;
    auto crvccch = event.getValidHandle<CrvCoincidenceClusterCollection>(_crvccc);
    auto const& crvccc = *crvccch;
    auto crvdch = event.getValidHandle<CrvDigiCollection>(_crvdc);
    auto const& crvdc = *crvdch;
    auto crvdmcch = event.getValidHandle<CrvDigiMCCollection>(_crvdmcc);
    auto const& crvdmcc = *crvdmcch;
    auto vdspch = event.getValidHandle<StepPointMCCollection>(_vdspc);
    auto const& vdspc = *vdspch;
    auto cdch = event.getValidHandle<CaloDigiCollection>(_cdc);
    auto const& cdc = *cdch;
    auto ccch = event.getValidHandle<CaloClusterCollection>(_ccc);
    auto const& ccc = *ccch;
    auto cssch = event.getValidHandle<CaloShowerSimCollection>(_cssc);
    auto const& cssc = *cssch;
    // some things needed for creating Ptrs before the collection is in the event
    auto KalSeedMCCollectionPID = getProductID<KalSeedMCCollection>();
    auto KalSeedMCCollectionGetter = event.productGetter(KalSeedMCCollectionPID);
    auto CaloClusterMCCollectionPID = getProductID<CaloClusterMCCollection>();
    auto CaloClusterMCCollectionGetter = event.productGetter(CaloClusterMCCollectionPID);
// create output; these are pruned collections containing only
// products related to the reconstruction output or the event primary
    std::unique_ptr<KalSeedMCCollection> ksmcc(new KalSeedMCCollection);
    std::unique_ptr<KalSeedMCAssns> ksmca(new KalSeedMCAssns);
    std::unique_ptr<CaloClusterMCCollection> ccmcc(new CaloClusterMCCollection);
    std::unique_ptr<CaloClusterMCAssns> ccmca(new CaloClusterMCAssns);
    std::unique_ptr<StrawDigiCollection> ssdc(new StrawDigiCollection);
    std::unique_ptr<StrawHitFlagCollection> sshfc(new StrawHitFlagCollection);
    std::unique_ptr<CaloDigiCollection> scdc(new CaloDigiCollection);
    std::unique_ptr<CrvDigiCollection> scrvdc(new CrvDigiCollection);
// index maps between original collections and pruned collections
    std::unique_ptr<IndexMap> sdmcim(new IndexMap);
    std::unique_ptr<IndexMap> crvdmcim(new IndexMap);
    std::unique_ptr<IndexMap> cdim(new IndexMap);
    // reco summary
    std::unique_ptr<RecoCount> nrec(new RecoCount);
    // straw hit (indices) that are referenced by the tracks, or particles
    // that contributed to the track
    SHIS shindices;
    // set of CaloClusters referenced by the tracks
    std::set<art::Ptr<CaloCluster> > ccptrs;
    // set of SimParticles to save for this event
    std::set<art::Ptr<SimParticle> > simps;
    // add the MC primary SimParticles
    for(auto const& spp : pp.primarySimParticles()) simps.insert(spp);
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
	  if(spcc.size() > 0)fillVDSP(det,spcc.front()._spp,vdspc,mcseed);
	  ksmcc->push_back(mcseed);
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
    // fill the StrawIndex map with the complete list of indices.
    StrawHitIndex shcount(0);
    ssdc->reserve(shindices.size());
    for(auto shindex : shindices){
      sdmcim->addElement(shindex,shcount++);
// deep-copy the selected StrawDigis and StrawHitFlags
      ssdc->push_back(sdc[shindex]);
      sshfc->push_back(shfc[shindex]);
    }
    if(_debug > 1) std::cout << "Selected " << shcount << " StrawDigis" << std::endl;
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
    for(auto crvindex : crvindices){
      crvdmcim->addElement(crvindex,crvcount++);
      // deep-copy the selected CrvDigis
      scrvdc->push_back(crvdc.at(crvindex));
    }
    // fill CaloClusters
    for(auto const& ccptr : ccptrs) {
      CaloClusterMC ccmc;
      fillCaloClusterMC(*ccptr,cssc,pp,ccmc);
      // save and create the assn
      ccmcc->push_back(ccmc);
      auto ccmcp = art::Ptr<CaloClusterMC>(CaloClusterMCCollectionPID,ccmcc->size()-1,CaloClusterMCCollectionGetter);
      ccmca->addSingle(ccptr,ccmcp);
  // deep-copy CaloDigis used in clusters
      for(auto const& cchptr : ccptr->caloCrystalHitsPtrVector()){
	for (auto const& rcdptr : cchptr->recoCaloDigis()){
	  scdc->push_back(*rcdptr->caloDigiPtr());
	}
      }
    }
    // Should add the SimParticles from the CaloClusters and CRV Coincidences FIXME!!

    // fill detailed StrawHit counts
    fillStrawHitCounts(shfc,*nrec.get());
    // general reco count
    nrec->_ncalodigi = cdc.size();
    nrec->_ncrvdigi = crvdc.size();
    nrec->_ncaloclust = ccc.size();
    // put output in event
    event.put(std::move(sdmcim),"StrawDigiMap");
    event.put(std::move(crvdmcim),"CrvDigiMap");
    event.put(std::move(cdim),"CaloDigiMap");
    event.put(std::move(ssdc));
    event.put(std::move(sshfc));
    event.put(std::move(scdc));
    event.put(std::move(scrvdc));
    event.put(std::move(ksmcc));
    event.put(std::move(ksmca));
    event.put(std::move(ccmcc));
    event.put(std::move(ccmca));
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
	if(sdmc.earlyStepPointMC()->simParticle() == spc._spp){
	  spref = isp;
	  break;
	}
      }
      if(spref < 0)throw cet::exception("Reco")<<"mu2e::SelectRecoMC: missing index"<< std::endl;
      tshmc._spindex = spref;
      // fill other info directly from the StrawDigiMC
      tshmc._energySum = sdmc.triggerEnergySum(sdmc.earlyEnd());
      const auto& mcstep = *(sdmc.earlyStepPointMC());
      tshmc._cpos = Geom::toXYZVec(sdmc.clusterPosition(sdmc.earlyEnd()));
      tshmc._mom = mcstep.momentum();
      tshmc._time = fmod(_toff.timeWithOffsetsApplied(mcstep),_mbtime);
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
	if(sdmc.earlyStepPointMC()->simParticle() == spc._spp){
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
	    const auto& mcstep = *(sdmc.earlyStepPointMC());
	    tshmc._cpos = Geom::toXYZVec(sdmc.clusterPosition(sdmc.earlyEnd()));
	    tshmc._mom = mcstep.momentum();
	    tshmc._time = fmod(_toff.timeWithOffsetsApplied(mcstep),_mbtime);
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

  void SelectRecoMC::fillCaloClusterMC(CaloCluster const& cc,
      CaloShowerSimCollection const& cssc,PrimaryParticle const& pp, 
      CaloClusterMC& ccmc){
    // struct for sorting MC energy deposits, largest first
    struct esort : public std::binary_function <CaloMCEDep, CaloMCEDep, bool> {
      bool operator() (CaloMCEDep a, CaloMCEDep b) { return a._edep > b._edep; }
    };
    //  keep track of energy contributions
    std::vector<CaloMCEDep> edeps;
    // loop over the crystal hits in this cluster
    for(auto const& cchptr : cc.caloCrystalHitsPtrVector()){
      // look for matches in the CaloShowerSimCollection: this is approximate, as the fundamental connectoion back
      // to CaloDigis is missing FIXME!
      // Must truncate tiny energy deposts as the CaloShowerSim compression sometimes fails!
      for(auto const& css : cssc) {
	if(css.energyMC() > _ccmce && css.crystalId() == cchptr->id()){
	  // compare times
	  double csstime = css.time();
	  if(_debug > 2) std::cout << "Matching CaloShowerSim crystal id " << css.crystalId()
	    << " MCenergy " << css.energyMC() << " time " << csstime 
	      << " calohit time " << cchptr->time() << std::endl;
	  if(fabs(csstime-cchptr->time()) < _ccmcdt) {
	    // match: see if an entry for this SimParticle already exists
	    bool found(false);
	    for(auto& edep : edeps) {
	      if(css.sim() == edep.sim()){
		// add in the energy; average the time
		float te = edep.energyDeposit()*edep.time() + csstime*css.energyMC();
		edep._edep += css.energyMC();
		edep._time = te/edep.energyDeposit();
		found = true;
		break;
	      }
	    }
	    if(!found){
	      // add a new element
	      CaloMCEDep edep;
	      edep._simp = css.sim();
	      edep._edep = css.energyMC();
	      edep._time = csstime;
	      // try all primary particle SimParticles, and keep the one with the best
	      // relationship
	      for(auto const& spp : pp.primarySimParticles()){
		MCRelationship mcr(spp,css.sim());
		if(mcr > edep._rel)edep._rel = mcr;
	      }
	      edeps.push_back(edep);
	    }
	  }
	}
      }
    }
    // sort
    std::sort(edeps.begin(),edeps.end(),esort());
    // fill the CaloClusterMC object
    ccmc._edep = 0.0;
    ccmc._time = 0.0;
    ccmc._edeps.reserve(edeps.size());
    for(auto& edep : edeps) {
      ccmc._edeps.push_back(edep);
      ccmc._edep += edep.energyDeposit();
      ccmc._time += edep.time()*edep.energyDeposit();
    }
    ccmc._time /= ccmc._edep;
  }

  void SelectRecoMC::fillStrawHitCounts(StrawHitFlagCollection const& shfc, RecoCount& nrec) {
    nrec._nstrawdigi = shfc.size();
    for(const auto& shf : shfc) {
      if(shf.hasAllProperties(StrawHitFlag::energysel))++nrec._nshfesel;
      if(shf.hasAllProperties(StrawHitFlag::radsel))++nrec._nshfrsel;
      if(shf.hasAllProperties(StrawHitFlag::timesel))++nrec._nshftsel;
      if(shf.hasAllProperties(StrawHitFlag::bkg))++nrec._nshfbkg;
      if(shf.hasAllProperties(StrawHitFlag::trksel))++nrec._nshftpk;
    }
  }
}
DEFINE_ART_MODULE(mu2e::SelectRecoMC)
