//
// Prototype analysis module using tracks.  This module associates information from the
// Mu2e detector systems into a single coherent analysis TTree (trkana).  This module
// depends on the data products produced by reconstruction and (if requested) their MC
// counterparts.  The primary analysis object is the set of Downstream electron track fits.
// Upstream electron fit and downstream muon are also required for PID and quality selection.
// Calorimeter clusters and Track-cluster matching are used for PID. CRV coincidences are also
// included for rejecting cosmic backgrounds.
// Most of the calcluations are done by upstream modules and helper classes.
// Original author: Dave Brown (LBNL) 7/7/2016
// Updated November 2018 to run on KalSeeds only (A. Edmonds)
//

// Mu2e includes
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/KalSeedAssns.hh"
#include "RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "GeometryService/inc/VirtualDetector.hh"
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"

// ROOT incldues
#include "Rtypes.h"
#include "TBits.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1F.h"

// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// mu2e tracking
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
// diagnostics
#include "TrkDiag/inc/TrkComp.hh"
#include "TrkDiag/inc/HitCount.hh"
#include "TrkDiag/inc/TrkCount.hh"
#include "TrkDiag/inc/EventInfo.hh"
#include "TrkDiag/inc/TrkInfo.hh"
#include "TrkDiag/inc/GenInfo.hh"
#include "TrkDiag/inc/EventWeightInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfo.hh"
#include "TrkDiag/inc/TrkStrawHitInfoMC.hh"
#include "TrkDiag/inc/TrkCaloHitInfo.hh"
#include "TrkDiag/inc/CaloClusterInfoMC.hh"
#include "TrkDiag/inc/TrkQualInfo.hh"
#include "TrkDiag/inc/TrkPIDInfo.hh"
#include "TrkDiag/inc/HelixInfo.hh"
#include "TrkDiag/inc/InfoStructHelper.hh"
#include "TrkDiag/inc/InfoMCStructHelper.hh"
#include "RecoDataProducts/inc/RecoQual.hh"
#include "TrkDiag/inc/RecoQualInfo.hh"
// CRV info
#include "CRVAnalysis/inc/CRVAnalysis.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

namespace mu2e {
// Need this for the BaBar headers.
  using CLHEP::Hep3Vector;
  typedef KalSeedCollection::const_iterator KSCIter;

  class TrackAnalysisReco : public art::EDAnalyzer {

  public:

    struct BranchOptConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<bool> fillmc{Name("fillMC"), Comment("Switch to turn on filling of MC information for this set of tracks"), false};
      fhicl::Atom<bool> fillhits{Name("fillHits"), Comment("Switch to turn on filling of hit-level information for this set of tracks"), false};
      fhicl::OptionalAtom<std::string> trkqual{Name("trkqual"), Comment("TrkQualCollection input tag to be written out (use prefix if fcl parameter suffix is defined)")};
      fhicl::Atom<bool> filltrkqual{Name("fillTrkQual"), Comment("Switch to turn on filling of the full TrkQualInfo for this set of tracks"), false};
      fhicl::OptionalAtom<std::string> trkpid{Name("trkpid"), Comment("TrkCaloHitPIDCollection input tag to be written out (use prefix if fcl parameter suffix is defined)")};
      fhicl::Atom<bool> filltrkpid{Name("fillTrkPID"), Comment("Switch to turn on filling of the full TrkPIDInfo for this set of tracks"), false};
      fhicl::Atom<bool> required{Name("required"), Comment("True/false if you require this type of track in the event"), false};
    };

    struct BranchConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<std::string> input{Name("input"), Comment("KalSeedCollection input tag (use prefix if fcl parameter suffix is defined)")};
      fhicl::Atom<std::string> branch{Name("branch"), Comment("Name of output branch")};
      fhicl::Atom<std::string> suffix{Name("suffix"), Comment("Fit suffix (e.g. DeM)"), ""};
      fhicl::Table<BranchOptConfig> options{Name("options"), Comment("Optional arguments for a branch")};
    };

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Table<BranchConfig> candidate{Name("candidate"), Comment("Candidate physics track info")};
      fhicl::OptionalSequence< fhicl::Table<BranchConfig> > supplements{Name("supplements"), Comment("Supplemental physics track info (TrkAna will find closest in time to candidate)")};
      fhicl::Atom<art::InputTag> rctag{Name("RecoCountTag"), Comment("RecoCount"), art::InputTag()};
      fhicl::Atom<art::InputTag> meanPBItag{Name("MeanBeamIntensity"), Comment("Tag for MeanBeamIntensity"), art::InputTag()};
      fhicl::Atom<art::InputTag> PBIwtTag{Name("PBIWeightTag"), Comment("Tag for PBIWeight") ,art::InputTag()};
      fhicl::Atom<art::InputTag> caloClusterMCTag{Name("CaloClusterMCTag"), Comment("Tag for CaloClusterMCCollection") ,art::InputTag()};
      fhicl::Atom<std::string> crvCoincidenceModuleLabel{Name("CrvCoincidenceModuleLabel"), Comment("CrvCoincidenceModuleLabel")};
      fhicl::Atom<std::string> crvCoincidenceMCModuleLabel{Name("CrvCoincidenceMCModuleLabel"), Comment("CrvCoincidenceMCModuleLabel")};
      fhicl::Atom<std::string> crvRecoPulseLabel{Name("CrvRecoPulseLabel"), Comment("CrvRecoPulseLabel")};
      fhicl::Atom<std::string> crvStepLabel{Name("CrvStepLabel"), Comment("CrvStepLabel")};
      fhicl::Atom<std::string> simParticleLabel{Name("SimParticleLabel"), Comment("SimParticleLabel")};
      fhicl::Atom<std::string> mcTrajectoryLabel{Name("MCTrajectoryLabel"), Comment("MCTrajectoryLabel")};
      fhicl::Atom<double> crvPlaneY{Name("CrvPlaneY"),2751.485};  //y of center of the top layer of the CRV-T counters
      fhicl::Table<SimParticleTimeOffset::Config> timeOffsets{ Name("TimeOffsets"), Comment("Time maps") };
      fhicl::Atom<std::string> crvWaveformsModuleLabel{ Name("CrvWaveformsModuleLabel"), Comment("CrvWaveformsModuleLabel")};
      fhicl::Atom<std::string> crvDigiModuleLabel{ Name("CrvDigiModuleLabel"), Comment("CrvDigiModuleLabel")};
      fhicl::Atom<bool> fillmc{Name("FillMCInfo"),Comment("Global switch to turn on/off MC info"),true};
      fhicl::Atom<bool> pempty{Name("ProcessEmptyEvents"),false};
      fhicl::Atom<bool> crv{Name("AnalyzeCRV"),false};
      fhicl::Atom<bool> crvpulses{Name("AnalyzeCRVPulses"),false};
      fhicl::Atom<bool> helices{Name("FillHelixInfo"),false};
      fhicl::Atom<bool> filltrkqual{Name("FillTrkQualInfo"),false};
      fhicl::Atom<bool> filltrkpid{Name("FillTrkPIDInfo"),false};
      fhicl::Atom<bool> filltrig{Name("FillTriggerInfo"),false};
      fhicl::Atom<std::string> trigpathsuffix{Name("TriggerPathSuffix"), "_trigger"}; // all trigger paths have this in the name
      fhicl::Atom<int> diag{Name("diagLevel"),1};
      fhicl::Atom<bool> fillhits{Name("FillHitInfo"),Comment("Global switch to turn on/off hit-level info"), false};
      fhicl::Atom<int> debug{Name("debugLevel"),0};
      fhicl::Atom<art::InputTag> primaryParticleTag{Name("PrimaryParticleTag"), Comment("Tag for PrimaryParticle"), art::InputTag()};
      fhicl::Atom<art::InputTag> kalSeedMCTag{Name("KalSeedMCAssns"), Comment("Tag for KalSeedMCAssn"), art::InputTag()};
      fhicl::Table<InfoMCStructHelper::Config> infoMCStructHelper{Name("InfoMCStructHelper"), Comment("Configuration for the InfoMCStructHelper")};
      fhicl::Atom<bool> fillmcxtra{Name("FillExtraMCSteps"),false};
      fhicl::OptionalSequence<art::InputTag> mcxtratags{Name("ExtraMCStepCollectionTags"), Comment("Input tags for any other StepPointMCCollections you want written out")};
      fhicl::OptionalSequence<std::string> mcxtrasuffix{Name("ExtraMCStepBranchSuffix"), Comment("The suffix to the branch for the extra MC steps (e.g. putting \"ipa\" will give a branch \"demcipa\")")};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit TrackAnalysisReco(const Parameters& conf);
    virtual ~TrackAnalysisReco() { }

    void beginJob();
    void beginSubRun(const art::SubRun & subrun ) override;
    void analyze(const art::Event& e);

  private:

    Config _conf;
    std::vector<BranchConfig> _allBranches; // candidates + supplements
    size_t _candidateIndex; // location in above vector that contains the candidate

    // track comparator
    TrkComp _tcomp;

    // main TTree
    TTree* _trkana;
    TProfile* _tht; // profile plot of track hit times: just an example
    // general event info branch
    double _meanPBI;
    EventInfo _einfo;
    // hit counting
    HitCount _hcnt;
    // track counting
    TrkCount _tcnt;
    // track branches (inputs)
    std::vector<art::Handle<KalSeedCollection> > _allKSCHs;
    // track branches (outputs)
    std::vector<TrkInfo> _allTIs;
    std::vector<TrkFitInfo> _allEntTIs, _allMidTIs, _allXitTIs;
    std::vector<TrkCaloHitInfo> _allTCHIs;
    // quality branches (inputs)
    std::vector<std::vector<art::Handle<RecoQualCollection> > > _allRQCHs; // outer vector is for each candidate/supplement, inner vector is all RecoQuals
    std::vector<art::Handle<TrkQualCollection> > _allTQCHs; // we will only allow one TrkQual object per candidate/supplement to be fully written out
    std::vector<art::Handle<TrkCaloHitPIDCollection> > _allTCHPCHs; // we will only allow one TrkCaloHitPID object per candidate/supplement to be fully written out
    // quality branches (outputs)
    std::vector<RecoQualInfo> _allRQIs;
    std::vector<TrkQualInfo> _allTQIs;
    std::vector<TrkPIDInfo> _allTPIs;
    // trigger information
    unsigned _trigbits;
    TH1F* _trigbitsh; // plot of trigger bits: just an example
    std::map<size_t,unsigned> _tmap; // map between path and trigger ID.  ID should come from trigger itself FIXME!
    // MC truth branches (inputs)
    art::Handle<PrimaryParticle> _pph;
    art::Handle<KalSeedMCAssns> _ksmcah;
    art::Handle<CaloClusterMCCollection> _ccmcch;
    art::InputTag _primaryParticleTag;
    std::vector<int> _entvids, _midvids, _xitvids;
    // MC truth branches (outputs)
    std::vector<TrkInfoMC> _allMCTIs;
    std::vector<GenInfo> _allMCGenTIs, _allMCPriTIs;
    std::vector<TrkInfoMCStep> _allMCEntTIs, _allMCMidTIs, _allMCXitTIs;
    std::vector<CaloClusterInfoMC> _allMCTCHIs;

    // hit level info branches 
    std::vector<std::vector<TrkStrawHitInfo>> _allTSHIs;
    std::vector<std::vector<TrkStrawMatInfo>> _allTSMIs;
    std::vector<std::vector<TrkStrawHitInfoMC>> _allTSHIMCs;

    // event weights
    std::vector<art::Handle<EventWeight> > _wtHandles;
    EventWeightInfo _wtinfo;
    // CRV info
    std::vector<CrvHitInfoReco> _crvinfo;
    int _bestcrv;
    std::vector<CrvHitInfoMC> _crvinfomc;
    CrvSummaryReco _crvsummary;
    CrvSummaryMC   _crvsummarymc;
    std::vector<CrvPlaneInfoMC> _crvinfomcplane;
    std::vector<CrvPulseInfoReco> _crvpulseinfo;
    std::vector<CrvWaveformInfo> _crvwaveforminfo;
    std::vector<CrvHitInfoMC> _crvpulseinfomc;
    // helices
    HelixInfo _hinfo;
    // struct helpers
    InfoStructHelper _infoStructHelper;
    InfoMCStructHelper _infoMCStructHelper;

    // helper functions
    void fillEventInfo(const art::Event& event);
    void fillTriggerBits(const art::Event& event,std::string const& process);
    void resetBranches();
    size_t findSupplementTrack(KalSeedCollection const& kcol,KalSeed const& candidate, bool sameColl);
    void fillAllInfos(const art::Handle<KalSeedCollection>& ksch, size_t i_branch, size_t i_kseed);

    template <typename T, typename TI>
    std::vector<art::Handle<T> > createSpecialBranch(const art::Event& event, const std::string& branchname, 
						     std::vector<art::Handle<T> >& handles, TI& infostruct, const std::string& selection = "");

};

  TrackAnalysisReco::TrackAnalysisReco(const Parameters& conf):
    art::EDAnalyzer(conf),
    _conf(conf()),
    _trigbitsh(0),
    _infoMCStructHelper(conf().infoMCStructHelper())
  {
    _midvids.push_back(VirtualDetectorId::TT_Mid);
    _midvids.push_back(VirtualDetectorId::TT_MidInner);
    _entvids.push_back(VirtualDetectorId::TT_FrontHollow);
    _entvids.push_back(VirtualDetectorId::TT_FrontPA);
    _xitvids.push_back(VirtualDetectorId::TT_Back);

    // collect both candidate and supplement branches into one place
    _allBranches.push_back(_conf.candidate());
    _candidateIndex = 0;
    std::vector<BranchConfig> supps;
    if (_conf.supplements(supps)) {
      for(const auto& i_supp : supps) {
	_allBranches.push_back(i_supp);
      }
    }

    // Create all the info structs
    for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
      TrkInfo ti;
      _allTIs.push_back(ti);
      TrkFitInfo ent, mid, xit;
      _allEntTIs.push_back(ent);
      _allMidTIs.push_back(mid);
      _allXitTIs.push_back(xit);
      
      TrkCaloHitInfo tchi;
      _allTCHIs.push_back(tchi);
      
      TrkInfoMC mcti;
      _allMCTIs.push_back(mcti);
      GenInfo mcgen, mcpri;
      _allMCGenTIs.push_back(mcgen);
      _allMCPriTIs.push_back(mcpri);
      TrkInfoMCStep mcent, mcmid, mcxit;
      _allMCEntTIs.push_back(mcent);
      _allMCMidTIs.push_back(mcmid);
      _allMCXitTIs.push_back(mcxit);
      
      CaloClusterInfoMC mctchi;
      _allMCTCHIs.push_back(mctchi);

      RecoQualInfo rqi;
      _allRQIs.push_back(rqi);
      TrkQualInfo tqi;
      _allTQIs.push_back(tqi);
      TrkPIDInfo tpi;
      _allTPIs.push_back(tpi);

      std::vector<TrkStrawHitInfo> tshi;
      _allTSHIs.push_back(tshi);
      std::vector<TrkStrawMatInfo> tsmi;
      _allTSMIs.push_back(tsmi);
      std::vector<TrkStrawHitInfoMC> tshimc;
      _allTSHIMCs.push_back(tshimc);
    }
  }

  void TrackAnalysisReco::beginJob( ){
    art::ServiceHandle<art::TFileService> tfs;
// create TTree
    _trkana=tfs->make<TTree>("trkana","track analysis");
    _tht=tfs->make<TProfile>("tht","Track Hit Time Profile",RecoCount::_nshtbins,-25.0,1725.0);
// add event info branch
    _trkana->Branch("evtinfo.",&_einfo,EventInfo::leafnames().c_str());
// hit counting branch
    _trkana->Branch("hcnt.",&_hcnt,HitCount::leafnames().c_str());
// track counting branch
    std::vector<std::string> trkcntleaves;
    for (const auto& i_branchConfig : _allBranches) {
      trkcntleaves.push_back(i_branchConfig.branch());
    }
    _trkana->Branch("tcnt",&_tcnt,_tcnt.leafnames(trkcntleaves).c_str());

// create all candidate and supplement branches
    for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
      BranchConfig i_branchConfig = _allBranches.at(i_branch);
      std::string branch = i_branchConfig.branch();
      _trkana->Branch(branch.c_str(),&_allTIs.at(i_branch),TrkInfo::leafnames().c_str());
      _trkana->Branch((branch+"ent").c_str(),&_allEntTIs.at(i_branch),TrkFitInfo::leafnames().c_str());
      _trkana->Branch((branch+"mid").c_str(),&_allMidTIs.at(i_branch),TrkFitInfo::leafnames().c_str());
      _trkana->Branch((branch+"xit").c_str(),&_allXitTIs.at(i_branch),TrkFitInfo::leafnames().c_str());
      _trkana->Branch((branch+"tch").c_str(),&_allTCHIs.at(i_branch),TrkCaloHitInfo::leafnames().c_str());
      if (_conf.filltrkqual() && i_branchConfig.options().filltrkqual()) {
	_trkana->Branch((branch+"trkqual").c_str(), &_allTQIs.at(i_branch), TrkQualInfo::leafnames().c_str());
      }
      if (_conf.filltrkpid() && i_branchConfig.options().filltrkpid()) {
	_trkana->Branch((branch+"trkpid").c_str(), &_allTPIs.at(i_branch), TrkPIDInfo::leafnames().c_str());
      }
      // optionally add hit-level branches
      // (for the time being diagLevel : 2 will still work, but I propose removing this at some point)
      if(_conf.diag() > 1 || (_conf.fillhits() && i_branchConfig.options().fillhits())){ 
	_trkana->Branch((branch+"tsh").c_str(),&_allTSHIs.at(i_branch));
	_trkana->Branch((branch+"tsm").c_str(),&_allTSMIs.at(i_branch));
      }
      // optionall add MC branches
      if(_conf.fillmc() && i_branchConfig.options().fillmc()){
	_trkana->Branch((branch+"mc").c_str(),&_allMCTIs.at(i_branch),TrkInfoMC::leafnames().c_str());
	_trkana->Branch((branch+"mcgen").c_str(),&_allMCGenTIs.at(i_branch),GenInfo::leafnames().c_str());
	_trkana->Branch((branch+"mcpri").c_str(),&_allMCPriTIs.at(i_branch),GenInfo::leafnames().c_str());
	_trkana->Branch((branch+"mcent").c_str(),&_allMCEntTIs.at(i_branch),TrkInfoMCStep::leafnames().c_str());
	_trkana->Branch((branch+"mcmid").c_str(),&_allMCMidTIs.at(i_branch),TrkInfoMCStep::leafnames().c_str());
	_trkana->Branch((branch+"mcxit").c_str(),&_allMCXitTIs.at(i_branch),TrkInfoMCStep::leafnames().c_str());
	_trkana->Branch((branch+"tchmc").c_str(),&_allMCTCHIs.at(i_branch),CaloClusterInfoMC::leafnames().c_str());
	// at hit-level MC information
	// (for the time being diagLevel will still work, but I propose removing this at some point)
	if(_conf.diag() > 1 || (_conf.fillhits() && i_branchConfig.options().fillhits())){ 
	  _trkana->Branch((branch+"tshmc").c_str(),&_allTSHIMCs.at(i_branch));
	}
      }
    }
// trigger info.  Actual names should come from the BeginRun object FIXME
    if(_conf.filltrig()) {
      _trkana->Branch("trigbits",&_trigbits,"trigbits/i");
    }
// calorimeter information for the downstream electron track
// CRV info
    if(_conf.crv()) {
      _trkana->Branch("crvinfo",&_crvinfo);
      _trkana->Branch("crvsummary",&_crvsummary);
      _trkana->Branch("bestcrv",&_bestcrv,"bestcrv/I");
      if(_conf.crvpulses()) {
        _trkana->Branch("crvpulseinfo",&_crvpulseinfo);
        _trkana->Branch("crvwaveforminfo",&_crvwaveforminfo);
      }
      if(_conf.fillmc()){
	if(_conf.crv())
        {
          _trkana->Branch("crvinfomc",&_crvinfomc);
          _trkana->Branch("crvsummarymc",&_crvsummarymc);
          _trkana->Branch("crvinfomcplane",&_crvinfomcplane);
          if(_conf.crvpulses())
            _trkana->Branch("crvpulseinfomc",&_crvpulseinfomc);
        }
      }
    }
// helix info
   if(_conf.helices()) _trkana->Branch("helixinfo",&_hinfo,HelixInfo::leafnames().c_str());
  }

  void TrackAnalysisReco::beginSubRun(const art::SubRun & subrun ) {
    // mean number of protons on target
    art::Handle<ProtonBunchIntensity> PBIHandle;
    subrun.getByLabel(_conf.meanPBItag(), PBIHandle);
    if(PBIHandle.isValid())
      _meanPBI = PBIHandle->intensity();
    // get bfield
    _infoStructHelper.updateSubRun();
  }

  void TrackAnalysisReco::analyze(const art::Event& event) {
    // update timing maps for MC
    if(_conf.fillmc()){
      _infoMCStructHelper.updateEvent(event);
    }

    // need to create and define the event weight branch here because we only now know the EventWeight creating modules that have been run through the Event
    std::vector<art::Handle<EventWeight> > eventWeightHandles;
    _wtHandles = createSpecialBranch(event, "evtwt", eventWeightHandles, _wtinfo);

    std::string process = "Digitize"; // Digitization process is where the trigger is run
  /// Get the KalSeedCollections for both the candidate and all supplements
    _allKSCHs.clear();
    _allRQCHs.clear();
    _allTQCHs.clear();
    _allTCHPCHs.clear();
    // get Helix Assns
     art::Handle<KalHelixAssns> khaH;
    if(_conf.helices()){ // find associated Helices
      BranchConfig i_branchConfig = _allBranches.at(0);
      art::InputTag kalSeedInputTag = i_branchConfig.input() + i_branchConfig.suffix();
      event.getByLabel(kalSeedInputTag,khaH);
    }
    for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
      BranchConfig i_branchConfig = _allBranches.at(i_branch);
      art::Handle<KalSeedCollection> kalSeedCollHandle;
      art::InputTag kalSeedInputTag = i_branchConfig.input() + i_branchConfig.suffix();
      event.getByLabel(kalSeedInputTag,kalSeedCollHandle);
      _allKSCHs.push_back(kalSeedCollHandle);
      // also create the reco qual branches
      std::vector<art::Handle<RecoQualCollection> > recoQualCollHandles;
      std::vector<art::Handle<RecoQualCollection> > selectedRQCHs;
      selectedRQCHs = createSpecialBranch(event, i_branchConfig.branch()+"qual", recoQualCollHandles, _allRQIs.at(i_branch), i_branchConfig.suffix());
      for (const auto& i_selectedRQCH : selectedRQCHs) {
	if (i_selectedRQCH->size() != kalSeedCollHandle->size()) {
	  throw cet::exception("TrkAna") << "Sizes of KalSeedCollection and this RecoQualCollection are inconsistent (" << kalSeedCollHandle->size() << " and " << i_selectedRQCH->size() << " respectively)";
	}
      }
      _allRQCHs.push_back(selectedRQCHs);

      // TrkQual
      std::string i_trkqual_tag;
      art::Handle<TrkQualCollection> trkQualCollHandle;
      if (i_branchConfig.options().trkqual(i_trkqual_tag) && i_branchConfig.options().filltrkqual() && _conf.filltrkqual()) {
	art::InputTag trkQualInputTag = i_trkqual_tag + i_branchConfig.suffix();
	event.getByLabel(trkQualInputTag,trkQualCollHandle);
	if (trkQualCollHandle->size() != kalSeedCollHandle->size()) {
	  throw cet::exception("TrkAna") << "Sizes of KalSeedCollection and TrkQualCollection are inconsistent (" << kalSeedCollHandle->size() << " and " << trkQualCollHandle->size() << " respectively)";
	}
      }
      _allTQCHs.push_back(trkQualCollHandle);

      // TrkCaloHitPID
      std::string i_trkpid_tag;
      art::Handle<TrkCaloHitPIDCollection> trkpidCollHandle;
      if (i_branchConfig.options().trkpid(i_trkpid_tag) && i_branchConfig.options().filltrkpid() && _conf.filltrkpid()) {
	art::InputTag trkpidInputTag = i_trkpid_tag + i_branchConfig.suffix();
	event.getByLabel(trkpidInputTag,trkpidCollHandle);
	if (trkpidCollHandle->size() != kalSeedCollHandle->size()) {
	  throw cet::exception("TrkAna") << "Sizes of KalSeedCollection and TrkCaloHitPIDCollection are inconsistent (" << kalSeedCollHandle->size() << " and " << trkpidCollHandle->size() << " respectively)";
	}
      }
      _allTCHPCHs.push_back(trkpidCollHandle);
    }

    // general reco counts
    auto rch = event.getValidHandle<RecoCount>(_conf.rctag());
    auto const& rc = *rch;
    for(size_t ibin=0;ibin < rc._nshtbins; ++ibin){
      float time = rc._shthist.binMid(ibin);
      float count  = rc._shthist.binContents(ibin);
      _tht->Fill(time,count);
    }

    // trigger information
    if(_conf.filltrig()){
      fillTriggerBits(event,process);
    }
    // MC data
    if(_conf.fillmc()) { // get MC product collections
      event.getByLabel(_conf.primaryParticleTag(),_pph);
      event.getByLabel(_conf.kalSeedMCTag(),_ksmcah);
      event.getByLabel(_conf.caloClusterMCTag(),_ccmcch);
    }
    // reset event level structs
    _einfo.reset();
    _hcnt.reset();
    _tcnt.reset();
    _hinfo.reset();
    _wtinfo.reset();
    // reset
    resetBranches();
    // fill track counts
    for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
      _tcnt._counts[i_branch] = (_allKSCHs.at(i_branch))->size();
    }

    // fill event level info
    fillEventInfo(event);
    _infoStructHelper.fillHitCount(rc, _hcnt);

    // loop through all candidate tracks
    const auto& candidateKSCH = _allKSCHs.at(_candidateIndex);
    const auto& candidateKSC = *candidateKSCH;
    for (size_t i_kseed = 0; i_kseed < candidateKSC.size(); ++i_kseed) {

      bool skip_kseed = false; // there may be a reason we don't want to write this KalSeed out

      auto const& candidateKS = candidateKSC.at(i_kseed);
      fillAllInfos(candidateKSCH, _candidateIndex, i_kseed); // fill the info structs for the candidate
      if(_conf.helices()){
	auto const& khassns = khaH.product();
      // find the associated HelixSeed to this KalSeed using the assns.
	auto hptr = (*khassns)[i_kseed].second;
	_infoStructHelper.fillHelixInfo(hptr, _hinfo);
      }

      // Now loop through all the branches (both candidate + supplements)...
      for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
	if (i_branch == _candidateIndex) { // ...but actually ignore candidate
	  continue;
	}
	// check if supplement input collection is the same as the candidate input collections
	bool sameColl = false;
	if ( (_allBranches.at(_candidateIndex).input()+_allBranches.at(_candidateIndex).suffix()) 
	     == (_allBranches.at(i_branch).input()+_allBranches.at(i_branch).suffix()) ) {
	  sameColl = true;
	}
	const auto& i_supplementKSCH = _allKSCHs.at(i_branch);
	const auto& i_supplementKSC = *i_supplementKSCH;

	// If we require a supplement track of this type, and there are none...
	if (i_supplementKSC.size()==0 && _allBranches.at(i_branch).options().required()) {
	  skip_kseed = true; // ...skip this KalSeed
	}	  

	// find the supplement track closest in time
	auto i_supplementKS = findSupplementTrack(i_supplementKSC,candidateKS,sameColl);
	if(i_supplementKS < i_supplementKSC.size()) { 
	  fillAllInfos(_allKSCHs.at(i_branch), i_branch, i_supplementKS);
	}
      }
    
      if (skip_kseed) {
	continue;
      }

      // TODO we want MC information when we don't have a track
      // fill CRV info
      if(_conf.crv()){
	CRVAnalysis::FillCrvHitInfoCollections(_conf.crvCoincidenceModuleLabel(), _conf.crvCoincidenceMCModuleLabel(),
                                               _conf.crvRecoPulseLabel(), _conf.crvStepLabel(), _conf.simParticleLabel(), _conf.mcTrajectoryLabel(), event,
                                               _crvinfo, _crvinfomc, _crvsummary, _crvsummarymc, _crvinfomcplane, _conf.crvPlaneY());
        if(_conf.crvpulses())
          CRVAnalysis::FillCrvPulseInfoCollections(_conf.crvRecoPulseLabel(), _conf.crvWaveformsModuleLabel(), _conf.crvDigiModuleLabel(),
                                                   _infoMCStructHelper.getTimeMaps(), event, _crvpulseinfo, _crvpulseinfomc, _crvwaveforminfo);

//	find the best CRV match (closest in time)
	_bestcrv=-1;
	float mindt=1.0e9;
	float t0 = candidateKS.t0().t0();
	for(size_t icrv=0;icrv< _crvinfo.size(); ++icrv){
	  auto const& crvinfo = _crvinfo[icrv];
	  float dt = std::min(fabs(crvinfo._timeWindowStart-t0), fabs(crvinfo._timeWindowEnd-t0) );
	  if(dt < mindt){
	    mindt =dt;
	    _bestcrv = icrv;
	  }
	}
      }
      // fill this row in the TTree
      _trkana->Fill();
    }

    if(_conf.pempty() && candidateKSC.size()==0) { // if we want to process empty events
      _trkana->Fill();
    }
  }

  size_t TrackAnalysisReco::findSupplementTrack(KalSeedCollection const& kcol,const KalSeed& candidate, bool sameColl) {
    size_t retval = kcol.size();

    // loop over supplement tracks and find the closest
    double candidate_time = candidate.t0().t0();
    double closest_time = 999999999;
    for(auto i_kseed=kcol.begin(); i_kseed != kcol.end(); i_kseed++) {
      double supplement_time = i_kseed->t0().t0();
      if( fabs(supplement_time - candidate_time) < fabs(closest_time-candidate_time)) {
	if (sameColl && fabs(supplement_time - candidate_time)<1e-5) {
	  continue; // don't want the exact same track
	}
	closest_time = supplement_time;
	retval = i_kseed - kcol.begin();
      }
    }
    return retval;
  }

  void TrackAnalysisReco::fillEventInfo( const art::Event& event) {
    // fill basic event information
    _einfo._eventid = event.event();
    _einfo._runid = event.run();
    _einfo._subrunid = event.subRun();

    // get event weight products
    std::vector<Float_t> weights;
    for (const auto& i_weightHandle : _wtHandles) {
      double weight = i_weightHandle->weight();
      if (i_weightHandle.provenance()->moduleLabel() == _conf.PBIwtTag().label()) {
	if (_meanPBI > 0.0){
	  _einfo._nprotons = _meanPBI*weight;
	}
	else {
	  _einfo._nprotons = 1; // for non-background mixed jobs
	}
      }
      weights.push_back(weight);
    }
    _wtinfo.setWeights(weights);
  }

  void TrackAnalysisReco::fillTriggerBits(const art::Event& event,std::string const& process) {
    //get the TriggerResult from the process that created the KalFinalFit downstream collection
    art::InputTag const tag{Form("TriggerResults::%s", process.c_str())};
    auto trigResultsH = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults* trigResults = trigResultsH.product();
    TriggerResultsNavigator tnav(trigResults);
    _trigbits = 0;
   // setup the bin labels
    if(_trigbitsh == 0){ // is there a better way to do this?  I think not
      unsigned ntrig(0);
      unsigned npath = trigResults->size();
      for(size_t ipath=0;ipath < npath; ++ipath){
	if (tnav.getTrigPath(ipath).find(_conf.trigpathsuffix()) != std::string::npos) {
	  _tmap[ipath] = ntrig;
	  ntrig++;
	}
      }
      // build trigger histogram
      art::ServiceHandle<art::TFileService> tfs;
      _trigbitsh = tfs->make<TH1F>("trigbits","Trigger IDs",ntrig,-0.5,ntrig-0.5);
      for(size_t ipath=0;ipath < npath; ++ipath){
	auto ifnd = _tmap.find(ipath);
	if(ifnd != _tmap.end()){
	  _trigbitsh->GetXaxis()->SetBinLabel(ifnd->second+1,tnav.getTrigPath(ipath).c_str());
	}
      }
    }
    for(size_t ipath=0;ipath < trigResults->size(); ++ipath){
      if(trigResults->accept(ipath)) {
	auto ifnd = _tmap.find(ipath);
	if(ifnd != _tmap.end()){
	  unsigned itrig = ifnd->second;
	  _trigbitsh->Fill(itrig);
	  _trigbits |= 1 << itrig;
	  if(_conf.debug() > 1)
	    cout << "Trigger path " << tnav.getTrigPath(ipath) << " Trigger ID " << itrig << " returns " << trigResults->accept(ipath) << endl;
	}
      }
    }
    if(_conf.debug() > 0){
      cout << "Found TriggerResults for process " << process << " with " << trigResults->size() << " Lines"
	<< " trigger bits word " << _trigbits << endl;
      if(_conf.debug() > 1){
	TriggerResultsNavigator tnav(trigResults);
	tnav.print();
      }
    }
  }

  void TrackAnalysisReco::fillAllInfos(const art::Handle<KalSeedCollection>& ksch, size_t i_branch, size_t i_kseed) {

    const auto& kseed = ksch->at(i_kseed);
    BranchConfig branchConfig = _allBranches.at(i_branch);

    // get VD positions
    mu2e::GeomHandle<VirtualDetector> vdHandle;
    mu2e::GeomHandle<DetectorSystem> det;
    const XYZVec& entpos = XYZVec(det->toDetector(vdHandle->getGlobal(*_entvids.begin())));
    const XYZVec& midpos = XYZVec(det->toDetector(vdHandle->getGlobal(*_midvids.begin())));
    const XYZVec& xitpos = XYZVec(det->toDetector(vdHandle->getGlobal(*_xitvids.begin())));

    _infoStructHelper.fillTrkInfo(kseed,_allTIs.at(i_branch));
    _infoStructHelper.fillTrkFitInfo(kseed,_allEntTIs.at(i_branch),entpos);
    _infoStructHelper.fillTrkFitInfo(kseed,_allMidTIs.at(i_branch),midpos);
    _infoStructHelper.fillTrkFitInfo(kseed,_allXitTIs.at(i_branch),xitpos);
    //      _tcnt._overlaps[0] = _tcomp.nOverlap(kseed, kseed);

    if(_conf.diag() > 1 || (_conf.fillhits() && branchConfig.options().fillhits())){ // want hit level info
      _infoStructHelper.fillHitInfo(kseed, _allTSHIs.at(i_branch));
      _infoStructHelper.fillMatInfo(kseed, _allTSMIs.at(i_branch));
    }

// calorimeter info
    if (kseed.hasCaloCluster()) {
      _infoStructHelper.fillCaloHitInfo(kseed,  _allTCHIs.at(i_branch));
      _tcnt._ndec = 1; // only 1 possible calo hit at the moment
      // test
      if(_conf.debug()>0){
	auto const& tch = kseed.caloHit();
	auto const& cc = tch.caloCluster();
	std::cout << "CaloCluster has energy " << cc->energyDep()
		  << " +- " << cc->energyDepErr() << std::endl;
      }
    }

// all RecoQuals
    std::vector<Float_t> recoQuals; // for the output value
    for (const auto& i_recoQualHandle : _allRQCHs.at(i_branch)) {
      Float_t recoQual = i_recoQualHandle->at(i_kseed)._value;
      recoQuals.push_back(recoQual);
      Float_t recoQualCalib = i_recoQualHandle->at(i_kseed)._calib;
      recoQuals.push_back(recoQualCalib);
    }
    _allRQIs.at(i_branch).setQuals(recoQuals);
// TrkQual
    std::string trkqual_branch;
    if(_conf.filltrkqual() && branchConfig.options().filltrkqual() && branchConfig.options().trkqual(trkqual_branch)) { 
      const auto& trkQualCollHandle = _allTQCHs.at(i_branch);
      if (trkQualCollHandle.isValid()) { // we could have put an empty TrkQualCollection in, if we didn't want it
	const auto& trkQualColl = *trkQualCollHandle;
	const auto& trkQual = trkQualColl.at(i_kseed);
	_infoStructHelper.fillTrkQualInfo(trkQual, _allTQIs.at(i_branch));
      }
    }
// TrkCaloHitPID
    std::string trkpid_branch;
    if (_conf.filltrkpid() && branchConfig.options().filltrkpid() && branchConfig.options().trkpid(trkpid_branch)) {
      const auto& tchpcolH = _allTCHPCHs.at(i_branch);
      if (tchpcolH.isValid()) {
	const auto& tchpcol = *tchpcolH;
	auto const& tpid = tchpcol.at(i_kseed);
	_infoStructHelper.fillTrkPIDInfo(tpid, kseed, _allTPIs.at(i_branch));
      }
    }
// fill MC info associated with this track
    if(_conf.fillmc() && branchConfig.options().fillmc()) { 
      const PrimaryParticle& primary = *_pph;
      // use Assns interface to find the associated KalSeedMC; this uses ptrs
      auto kptr = art::Ptr<KalSeed>(ksch,i_kseed);
      //	std::cout << "KalSeedMCMatch has " << _ksmcah->size() << " entries" << std::endl;
      for(auto iksmca = _ksmcah->begin(); iksmca!= _ksmcah->end(); iksmca++){
	//	  std::cout << "KalSeed Ptr " << kptr << " match Ptr " << iksmca->first << std::endl;
	if(iksmca->first == kptr) {
	  auto const& kseedmc = *(iksmca->second);
	  _infoMCStructHelper.fillTrkInfoMC(kseedmc, _allMCTIs.at(i_branch));
	  double t0 = kseed.t0().t0();
	  _infoMCStructHelper.fillTrkInfoMCStep(kseedmc, _allMCEntTIs.at(i_branch), _entvids, t0);
	  _infoMCStructHelper.fillTrkInfoMCStep(kseedmc, _allMCMidTIs.at(i_branch), _midvids, t0);
	  _infoMCStructHelper.fillTrkInfoMCStep(kseedmc, _allMCXitTIs.at(i_branch), _xitvids, t0);
	  _infoMCStructHelper.fillGenAndPriInfo(kseedmc, primary, _allMCPriTIs.at(i_branch), _allMCGenTIs.at(i_branch));
	    
	  if(_conf.diag() > 1 || (_conf.fillhits() && branchConfig.options().fillhits())){ 
	    _infoMCStructHelper.fillHitInfoMCs(kseedmc, _allTSHIMCs.at(i_branch));
	  }
	  break;
	}
      }
      if (kseed.hasCaloCluster()) {
	// fill MC truth of the associated CaloCluster.  Use the fact that these are correlated by index with the clusters in that collection
	auto index = kseed.caloCluster().key();
	auto const& ccmcc = *_ccmcch;
	auto const& ccmc = ccmcc[index];
	_infoMCStructHelper.fillCaloClusterInfoMC(ccmc,_allMCTCHIs.at(i_branch));  // currently broken due to CaloMC changes.  This needs fixing in compression
      }
    }
  }

  // some branches can't be made until the analyze() function because we want to write out all data products of a certain type
  template <typename T, typename TI>
  std::vector<art::Handle<T> >  TrackAnalysisReco::createSpecialBranch(const art::Event& event, const std::string& branchname, 
								       std::vector<art::Handle<T> >& handles, TI& infostruct, const std::string& selection) {
    std::vector<art::Handle<T> > outputHandles = event.getMany<T>();
    if (handles.size()>0) {
      std::vector<std::string> labels;
      for (const auto& i_handle : handles) {
	std::string moduleLabel = i_handle.provenance()->moduleLabel();
	// event.getMany() doesn't have a way to wildcard part of the ModuleLabel, do it ourselves here
	size_t pos;
	if (selection != "") { // if we want to add a selection
	  pos = moduleLabel.find(selection);

	  // make sure that the selection (e.g. "DeM") appears at the end of the module label
	  if (pos == std::string::npos) {
	    //      std::cout << "Selection not found" << std::endl;
	    continue;
	  }
	  else if (pos+selection.length() != moduleLabel.size()) {
	    //      std::cout << "Selection wasn't at end of moduleLabel" << std::endl;
	    continue;
	  }
	  moduleLabel = moduleLabel.erase(pos, selection.length());
	}
	std::string instanceName = i_handle.provenance()->productInstanceName();

	std::string branchname = moduleLabel;
	if (instanceName != "") {
	  branchname += "_" + instanceName;
	}
	outputHandles.push_back(i_handle);
	labels.push_back(branchname);
      }
      if (!_trkana->GetBranch(branchname.c_str())) {  // only want to create the branch once
	_trkana->Branch(branchname.c_str(), &infostruct, infostruct.leafnames(labels).c_str());
      }
    }
    return outputHandles;
  }

  void TrackAnalysisReco::resetBranches() {
    for (size_t i_branch = 0; i_branch < _allBranches.size(); ++i_branch) {
      _allTIs.at(i_branch).reset();
      _allEntTIs.at(i_branch).reset();
      _allMidTIs.at(i_branch).reset();
      _allXitTIs.at(i_branch).reset();

      _allTCHIs.at(i_branch).reset();

      _allMCTIs.at(i_branch).reset();
      _allMCGenTIs.at(i_branch).reset();
      _allMCPriTIs.at(i_branch).reset();

      _allMCEntTIs.at(i_branch).reset();
      _allMCMidTIs.at(i_branch).reset();
      _allMCXitTIs.at(i_branch).reset();
      _allMCTCHIs.at(i_branch).reset();

      _allRQIs.at(i_branch).reset();
      _allTQIs.at(i_branch).reset();
      _allTPIs.at(i_branch).reset();

      // clear vectors
      _allTSHIs.at(i_branch).clear();
      _allTSMIs.at(i_branch).clear();
      _allTSHIMCs.at(i_branch).clear();
    }
// clear vectors
    _crvinfo.clear();
    _crvinfomc.clear();
    _crvinfomcplane.clear();
    _crvpulseinfo.clear();
    _crvwaveforminfo.clear();
    _crvpulseinfomc.clear();
  }
}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::TrackAnalysisReco;
DEFINE_ART_MODULE(TrackAnalysisReco);
