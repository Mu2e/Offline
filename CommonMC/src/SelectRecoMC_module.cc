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
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/DataProducts/inc/IndexMap.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/RecoCount.hh"
// Utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrkDiag/inc/TrkMCTools.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

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
        fhicl::Atom<bool> trkOnly                 { Name("TrkOnly"), Comment("Ignore cal and crv"), false};
        fhicl::Atom<bool> saveEnergySteps         { Name("SaveEnergySteps"),                  Comment("Save all StepPoints that contributed energy to a StrawDigi"), false};
        fhicl::Atom<bool> saveUnused              { Name("SaveUnusedDigiMCs"),                  Comment("Save StrawDigiMCs from particles used in any fit"), true};
        fhicl::Atom<bool> saveAllUnused           { Name("SaveAllUnusedDigiMCs"),                  Comment("Save all StrawDigiMCs from particles used in any fit"), false};
        fhicl::Atom<art::InputTag> PP             { Name("PrimaryParticle"),                  Comment("PrimaryParticle")};
        fhicl::Atom<art::InputTag> CCC            { Name("CaloClusterCollection"),          Comment("CaloClusterCollection")};
        fhicl::Sequence<art::InputTag> CrvCCCs    { Name("CrvCoincidenceClusterCollections"),Comment("CrvCoincidenceClusterCollections")};
        fhicl::Atom<art::InputTag> SDC            { Name("StrawDigiCollection"),                  Comment("StrawDigiCollection")};
        fhicl::Atom<art::InputTag> CHC            { Name("ComboHitCollection"),                  Comment("ComboHitCollection for the original StrawHits (not Panel hits)")};
        fhicl::Atom<art::InputTag> CDC            { Name("CaloDigiCollection"),                  Comment("CaloDigiCollection")};
        fhicl::Atom<art::InputTag> SDMCC          { Name("StrawDigiMCCollection"),          Comment("StrawDigiMCCollection")};
        fhicl::Atom<art::InputTag> CRVDC          { Name("CrvDigiCollection"),                  Comment("CrvDigiCollection")};
        fhicl::Atom<art::InputTag> CRVDMCC        { Name("CrvDigiMCCollection"),                  Comment("CrvDigiMCCollection")};
        fhicl::Atom<art::InputTag> PBTMC          { Name("PBTMC"),                  Comment("ProtonBunchTimeMC")};
        fhicl::Atom<art::InputTag> EWM            { Name("EventWindowMarker"), Comment("EventWindowMarker")};
        fhicl::Sequence<std::string> KalSeeds     { Name("KalSeedCollections"),                  Comment("KalSeedCollections")};
        fhicl::Sequence<std::string> HelixSeeds   { Name("HelixSeedCollections"),                  Comment("HelixSeedCollections")};
        fhicl::Atom<art::InputTag> VDSPC          { Name("VDSPCollection"),                  Comment("Virtual Detector StepPointMC collection")};
        fhicl::Atom<double> CCME                  { Name("CaloClusterMinE"),                Comment("Minimum energy CaloCluster to save digis (MeV)")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit SelectRecoMC(const Parameters& conf);
      void produce(art::Event& evt) override;

    private:
      typedef std::vector<TrkMCTools::spcount>SPCC;
      typedef std::set<StrawDigiIndex> SDIS;
      // utility functions
      void fillTSHMC         (KalSeed const& seed, SPCC const& spcc, StrawDigiMCCollection const& sdmcc, Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep, KalSeedMC& mcseed);
      void fillUnusedTSHMC   (SPCC const& spcc, StrawDigiMCCollection const& sdmcc, Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep, KalSeedMC& mcseed);
      void fillTSHMC        (TrkStrawHitMC& tshmc, size_t isdmc, size_t isp,StrawDigiMC const& sdmc, Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep );
      void fillVDSP          (GeomHandle<DetectorSystem>const& det, art::Ptr<SimParticle> const& psp, StepPointMCCollection const& vdspc, KalSeedMC& mcseed);
      void fillSPStubs       (SPCC const& spcc, PrimaryParticle const& pp, KalSeedMC& mcseed);
      void fillSDMCI         (KalSeedMC const& mcseed,SDIS& sdindices);
      void fillStrawHitCounts(ComboHitCollection const& chc, RecoCount& nrec);
      void fillTrk           (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, PrimaryParticle const& pp, RecoCount& nrec);
      void fillCrv           (art::Event& event, PrimaryParticle const& pp, RecoCount& nrec);
      void fillCalo          (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs,PrimaryParticle const& pp, RecoCount& nrec);
      int _debug;
      bool _trkonly;
      bool _saveallenergy, _saveunused, _saveallunused;
      art::InputTag _pp, _ccc, _sdc, _chc, _cdc, _sdmcc, _crvdc, _crvdmcc, _pbtmc, _ewm, _vdspc;
      std::vector<std::string> _kscs, _hscs;
    std::vector<art::InputTag> _crvcccs;
      double _ccme;
      // cache
      double _mbtime; // period of 1 microbunch
      double _pbtimemc; // mc true proton bunch time
      bool _onSpill;
      ProditionsHandle<StrawResponse> _strawResponse_h;
      ProditionsHandle<Tracker> _alignedTrackerSim_h{"Sim"};
  };

  SelectRecoMC::SelectRecoMC(const Parameters& config )  :
    art::EDProducer{config},
    _debug(config().debug()),
    _trkonly(config().trkOnly()),
    _saveallenergy(config().saveEnergySteps()),
    _saveunused(config().saveUnused()),
    _saveallunused(config().saveAllUnused()),
    _pp(config().PP()),
    _ccc(config().CCC()),
    _sdc(config().SDC()),
    _chc(config().CHC()),
    _cdc(config().CDC()),
    _sdmcc(config().SDMCC()),
    _crvdc(config().CRVDC()),
    _crvdmcc(config().CRVDMCC()),
    _pbtmc(config().PBTMC()),
    _ewm(config().EWM()),
    _vdspc(config().VDSPC()),
    _kscs(config().KalSeeds()),
    _hscs(config().HelixSeeds()),
    _crvcccs(config().CrvCCCs()),
    _ccme(config().CCME())
    {
      consumes<StrawDigiCollection>(_sdc);
      consumes<StrawDigiADCWaveformCollection>(_sdc);
      consumes<ComboHitCollection>(_chc);
      consumes<CaloDigiCollection>(_cdc);
      consumes<CaloClusterCollection>(_ccc);
      consumes<CrvDigiCollection>(_crvdc);
      consumesMany<KalSeedCollection>();
      for (const auto& i_tag : _crvcccs) {
        consumes<CrvCoincidenceClusterCollection>(i_tag);
      }
      consumes<PrimaryParticle>(_pp);
      consumes<StrawDigiMCCollection>(_sdmcc);
      consumes<CrvDigiMCCollection>(_crvdmcc);
      consumes<ProtonBunchTimeMC>(_pbtmc);
      consumes<EventWindowMarker>(_ewm);
      produces <IndexMap>("StrawDigiMap");
      if (!_trkonly){
        produces <IndexMap>("CrvDigiMap");
        produces <CaloDigiCollection>();
        produces <CrvDigiCollection>();
        produces <CrvRecoPulseCollection>();
        for (const auto& i_tag : _crvcccs) {
          produces <CrvCoincidenceClusterCollection>(i_tag.label());
        }
      }
      produces <KalSeedMCCollection>();
      produces <KalSeedMCAssns>();
      produces <StrawDigiCollection>();
      produces <StrawDigiADCWaveformCollection>();
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
    _mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
    auto pph = event.getValidHandle<PrimaryParticle>(_pp);
    auto const& pp = *pph;
    auto pbtmc = event.getValidHandle<ProtonBunchTimeMC>(_pbtmc);
    _pbtimemc = pbtmc->pbtime_;
    auto  ewmh = event.getValidHandle<EventWindowMarker>(_ewm);
    _onSpill = (ewmh->spillType() == EventWindowMarker::SpillType::onspill);

    std::unique_ptr<RecoCount> nrec(new RecoCount);
    std::set<art::Ptr<CaloCluster> > ccptrs;
    fillTrk(event,ccptrs,pp,*nrec.get());
    if (!_trkonly){
      fillCrv(event, pp, *nrec.get());
      fillCalo(event, ccptrs, pp, *nrec.get());
    }

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
      StrawDigiMCCollection const& sdmcc, Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep, KalSeedMC& mcseed) {
    for(auto const& hit : seed.hits() ) {
      // create a TrkStrawHitMC for each hit on the seed
      // find the referenced sim particle
      int spref(-1);
      auto const& sdmc = sdmcc.at(hit.index()); // bounds-check for security;
      // if mc info is not meaningful, do not try to inspect any SimParticles
      if (sdmc.provenance() != StrawDigiProvenance::External){
        for(size_t isp=0;isp < spcc.size(); isp++){
          auto const& spc = spcc[isp];
          if(sdmc.earlyStrawGasStep()->simParticle() == spc._spp){
            spref = isp;
            break;
          }
        }
        if(spref < 0)throw cet::exception("Reco")<<"mu2e::SelectRecoMC: missing index"<< std::endl;
      }
      TrkStrawHitMC tshmc;
      fillTSHMC(tshmc,hit.index(),spref,sdmc,tracker,srep);
      mcseed._tshmcs.push_back(tshmc);
    }
  }

  void SelectRecoMC::fillUnusedTSHMC( SPCC const& spcc, StrawDigiMCCollection const& sdmcc,
      Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep,
      KalSeedMC& mcseed) {
    // either keep hits only from the primary or all contributing particles,
    // assuming any exist. if they do not, then there is no primary.
    size_t limit = 0;
    if (0 < spcc.size()){
      limit = 1;
    }
    size_t ispmax = _saveallunused? spcc.size() : limit;
    for(size_t isp=0; isp < ispmax; ++isp){
      auto const& spc = spcc[isp];
      for (size_t isdmc=0; isdmc < sdmcc.size(); isdmc++){
        auto const& sdmc = sdmcc[isdmc];
        // if this contains no MC information, then we cannot inspect it further
        if (sdmc.provenance() == StrawDigiProvenance::External){
          continue;
        }
        auto const& sgs = *(sdmc.earlyStrawGasStep());
        if(sgs.simParticle() == spc._spp){
          // search to see if the associated digi is already on the track
          bool used(false);
          for(auto const& tshmc : mcseed._tshmcs ) {
            if(isdmc == tshmc.strawDigiMCIndex()){
              used = true;
              break;
            }
          }
          if(!used){
            TrkStrawHitMC tshmc;
            fillTSHMC(tshmc,isdmc,isp,sdmc,tracker,srep);
            mcseed._tshmcs.push_back(tshmc);
          }
        }
      }
    }
  }

  void SelectRecoMC::fillTSHMC(TrkStrawHitMC& tshmc, size_t isdmc, size_t isp,StrawDigiMC const& sdmc,
      Tracker const& tracker, std::shared_ptr<const StrawResponse>const& srep ) {
    // always record the index, to avoid downstream corruptions
    tshmc._sdmcindex = isdmc;
    // propagate interpretability of StrawDigiMC to TrkStrawHitMC
    if (sdmc.provenance() == StrawDigiProvenance::External){
      tshmc._provenance = TrkStrawHitProvenance::External;
      // if there is no MC information at all, leave TrkStrawHitMC empty
      return;
    }
    else if (sdmc.provenance() == StrawDigiProvenance::Mixed){
      tshmc._provenance = TrkStrawHitProvenance::Mixed;
    }
    else if (sdmc.provenance() == StrawDigiProvenance::Simulation){
      tshmc._provenance = TrkStrawHitProvenance::Simulation;
    }
    tshmc._spindex = isp;
    tshmc._energySum = sdmc.triggerEnergySum(sdmc.earlyEnd());
    const auto& sgs = *(sdmc.earlyStrawGasStep());
    tshmc._cpos = XYZVectorF(sdmc.clusterPosition(sdmc.earlyEnd()));
    tshmc._mom = sgs.momentum();
    tshmc._time = sgs.time();
    if (_onSpill){
      tshmc._time = fmod(tshmc._time,_mbtime);
      // fix for DAQ wrapping
      if(tshmc._time < -_pbtimemc)tshmc._time += _mbtime;
    }
    tshmc._strawId = sdmc.strawId();
    tshmc._earlyend = sdmc.earlyEnd();
    // compute the signal propagation time and drift time
    double vprop = 2.0*srep->halfPropV(sdmc.strawId(),1000.0*tshmc._energySum);
    auto const& straw = tracker.straw(sdmc.strawId());
    auto wdir = XYZVectorF(straw.wireDirection());
    auto tdir = sgs.momentum().Unit();
    double pdist = (straw.wireEnd(sdmc.earlyEnd())-sdmc.clusterPosition(sdmc.earlyEnd())).dot(straw.wireDirection());
    tshmc._tprop = fabs(pdist)/vprop;
    tshmc._tdrift = sdmc.wireEndTime(sdmc.earlyEnd()) -tshmc._time - tshmc._tprop - _pbtimemc - 2.4; // temporary kludge offset FIXME!
    if (_onSpill)
      tshmc._tdrift = fmod(tshmc._tdrift,_mbtime);
    const static XYZVectorF bdir(0.0,0.0,1.0);

    TwoLinePCA_XYZ wirepca( XYZVectorF(straw.wirePosition()), XYZVectorF(straw.wireDirection()),
          tshmc._cpos, tdir );
    TwoLinePCA_XYZ strawpca( XYZVectorF(straw.strawPosition()), XYZVectorF(straw.strawDirection()),
          tshmc._cpos, tdir );

    // vector from middle of wire to POCA
    auto wireMid_to_POCA = wirepca.point2() - XYZVectorF(straw.wirePosition());
    // longitudinal distance along wire to POCA
    tshmc._wireLen = wireMid_to_POCA.Dot(wdir);
    // perpedicular vector from wire to POCA
    auto POCA_delta = wirepca.point2() - wirepca.point1();
    // define DOCA to be negative in direction of track cross wire
    auto doca_sign_dir = tdir.Cross(wdir).Unit();
    tshmc._wireDOCA = -1*doca_sign_dir.Dot(POCA_delta);
    // phi is defined by vector to POCA, with 0 in 'Bdir'
    // and +pi/2 in V dir (radially out)
    auto Vdir = wdir.Cross(bdir);
    // ensure Vdir is pointing radially out
    if (Vdir.Dot(XYZVectorF(straw.wirePosition())) < 0.0) Vdir *= -1.0;
    tshmc._wirePhi = atan2(POCA_delta.Dot(Vdir),POCA_delta.Dot(bdir));
    // rdrift is the expected T2D given tdrift and phi
    tshmc._rdrift = srep->strawDrift().T2D(tshmc._tdrift,tshmc._wirePhi);
    // wireDot is cos angle between track and wire
    tshmc._wireDot = tdir.Dot(wdir);
    // wireTau is delta distance perpedicular to the wire from POCA to trigger cluster
    // first define unit vector perpendicular to particle track and wire direction
    auto track_cross_wire = tdir.Cross(wdir).Unit();
    // then this cross wdir is the vector perpendicular to the wire and the vector to POCA
    auto wire_cross_delta = wdir.Cross(track_cross_wire);
    // relative position of the trigger cluster
    auto wireMid_to_cluster = tshmc._cpos-XYZVectorF(straw.wirePosition());
    // wireTau is then the relative position of the trigger cluster in this direction
    tshmc._wireTau = wireMid_to_cluster.Dot(wire_cross_delta);

    auto sdir = XYZVectorF(straw.strawDirection());
    // perpedicular vector from straw center to POCA
    auto sPOCA_delta = strawpca.point2() - strawpca.point1();
    // define DOCA to be negative in direction of track cross straw
    auto sdoca_sign_dir = tdir.Cross(sdir).Unit();
    tshmc._strawDOCA = -1*sdoca_sign_dir.Dot(sPOCA_delta);
    // phi is defined by vector to POCA, with 0 in 'Bdir'
    // and +pi/2 in V dir (radially out)
    auto sVdir = sdir.Cross(bdir);
    // ensure Vdir is pointing radially out
    if (sVdir.Dot(XYZVectorF(straw.strawPosition())) < 0.0) sVdir *= -1.0;
    tshmc._strawPhi = atan2(sPOCA_delta.Dot(sVdir),sPOCA_delta.Dot(bdir));
  }

  void SelectRecoMC::fillSDMCI(KalSeedMC const& mcseed, SDIS& sdindices) {
    for(auto tshmc : mcseed.trkStrawHitMCs()) {
      sdindices.insert(tshmc.strawDigiMCIndex());
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
            vdsp.time(),
            vdsp.virtualDetectorId());
        mcseed._vdsteps.push_back(vds);
      }
    }
  }

  void SelectRecoMC::fillStrawHitCounts(ComboHitCollection const& chc, RecoCount& nrec) {

    for(const auto& ch : chc) {
      auto const& shf = ch.flag();
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
    auto srep = _strawResponse_h.getPtr(event.id());
    Tracker const& tracker = _alignedTrackerSim_h.get(event.id());

    // Tracker-reated data products
    auto sdch = event.getValidHandle<StrawDigiCollection>(_sdc);
    auto const& sdc = *sdch;
    auto sdadcch = event.getValidHandle<StrawDigiADCWaveformCollection>(_sdc);
    auto const& sdadcc = *sdadcch;
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
    // index maps between original collections and pruned collections
    std::unique_ptr<IndexMap> sdmcim(new IndexMap);
    // straw digi indices that are referenced by the tracks, or particles
    // that contributed to the track
    SDIS sdindices;
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
          fillTSHMC(seed,spcc,sdmcc,tracker,srep,mcseed);
          // add DigiMCs not used in the track but from true particles used in the track
          if(_saveunused)fillUnusedTSHMC(spcc,sdmcc,tracker,srep,mcseed);
          if(spcc.size() > 0 && vdspch.isValid()){
            auto const& vdspc = *vdspch;
            fillVDSP(det,spcc.front()._spp,vdspc,mcseed);
          }
          ksmcc->push_back(mcseed);
          // fill indices from all digis; those on the track and those from the MC true particle too
          fillSDMCI(mcseed,sdindices);
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
    // get digi indices from all helices too
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
          // go back to StrawDigi indices
          std::vector<StrawDigiIndex> shids;
          for(size_t ihit = 0; ihit < seed.hits().size(); ihit++)
            seed.hits().fillStrawDigiIndices(ihit,shids);
          // add these to the set (duplicates are suppressed)
          for(auto shid : shids)
            sdindices.insert(shid);
        }
      }
    }
    // fill the StrawIndex map with the complete list of indices.
    StrawDigiIndex sdcount(0);
    ssdc->reserve(sdindices.size());
    ssdadcc->reserve(sdindices.size());
    for(auto sdindex : sdindices){
      sdmcim->addElement(sdindex,sdcount++);
      // deep-copy the selected StrawDigis
      ssdc->push_back(sdc[sdindex]);
      ssdadcc->push_back(sdadcc[sdindex]);
    }
    if(_debug > 1) std::cout << "Selected " << sdcount << " StrawDigis" << std::endl;

    // fill detailed StrawHit counts
    fillStrawHitCounts(chc,nrec);
    event.put(std::move(sdmcim),"StrawDigiMap");
    event.put(std::move(ssdc));
    event.put(std::move(ssdadcc));
    event.put(std::move(ksmcc));
    event.put(std::move(ksmca));
  }

  void SelectRecoMC::fillCrv(art::Event& event,
      PrimaryParticle const& pp, RecoCount& nrec) {
    // find Crv data in event
    auto crvdch = event.getValidHandle<CrvDigiCollection>(_crvdc);
    auto const& crvdc = *crvdch;
    auto crvdmcch = event.getValidHandle<CrvDigiMCCollection>(_crvdmcc);
    auto const& crvdmcc = *crvdmcch;
    auto CrvRecoPulseCollectionPID = event.getProductID<CrvRecoPulseCollection>();
    auto CrvRecoPulseCollectionGetter = event.productGetter(CrvRecoPulseCollectionPID);
    // create new Crv collections
    std::unique_ptr<CrvDigiCollection> scrvdc(new CrvDigiCollection);
    std::unique_ptr<CrvRecoPulseCollection> scrvrpc(new CrvRecoPulseCollection);
    std::unique_ptr<IndexMap> crvdmcim(new IndexMap);
    std::vector<std::unique_ptr<CrvCoincidenceClusterCollection>> scrvcccs;

    std::set<uint16_t> crvindices;
    for (size_t i_tag = 0; i_tag < _crvcccs.size(); ++i_tag) {
      auto crvccch = event.getValidHandle<CrvCoincidenceClusterCollection>(_crvcccs.at(i_tag));
      auto const& crvccc = *crvccch;
      scrvcccs.push_back(std::unique_ptr<CrvCoincidenceClusterCollection>(new CrvCoincidenceClusterCollection));
      // loop over CrvCoincidenceClusters
      for(auto const& crvcc: crvccc) {
        std::vector<art::Ptr<CrvRecoPulse>> pulses;
        for(auto const& crvrp : crvcc.GetCrvRecoPulses()){
          // deep-copy the pulses used in coincidences: the digi indices are updated later
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
        scrvcccs.at(i_tag)->push_back(scrvcc);
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
    // update digi indices in the pulses
    for(auto& crvrp : *scrvrpc) {
      auto& indices = crvrp.GetWaveformIndices();
      for(size_t iindex=0; iindex < indices.size(); ++iindex){
        indices[iindex] = crvdmcim->getCondensedIndex(indices[iindex]);
      }
    }
    // update reco count
    nrec._ncrvdigi = crvdc.size();
    // put new data in event
    event.put(std::move(scrvdc));
    event.put(std::move(scrvrpc));
    for (size_t i_tag = 0; i_tag < _crvcccs.size(); ++i_tag) {
      event.put(std::move(scrvcccs.at(i_tag)), _crvcccs.at(i_tag).label());
    }
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
