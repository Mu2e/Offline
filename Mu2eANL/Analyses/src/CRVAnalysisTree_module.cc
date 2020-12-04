#include "TTree.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "TrkDiag/inc/EventInfo.hh"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "CRVAnalysis/inc/CRVAnalysis.hh"
#include "CLHEP/Random/Randomize.h"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"

namespace mu2e{
  class CRVAnalysisTree : public art::EDAnalyzer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<bool> crv{Name("AnalyzeCRV"),false};
      fhicl::Atom<bool> fillmc{Name("FillMCInfo"),true};
      fhicl::Atom<double> crvPlaneY{Name("CrvPlaneY"),2751.485};  //y of center of the top layer of the CRV-T counters
      fhicl::Atom<std::string> crvStepPointMCLabel{Name("CrvStepPointMCLabel"), Comment("CrvStepPointMCLabel")};
      fhicl::Atom<std::string> simParticleLabel{Name("SimParticleLabel"), Comment("SimParticleLabel")};
      fhicl::Atom<std::string> mcTrajectoryLabel{Name("MCTrajectoryLabel"), Comment("MCTrajectoryLabel")};

      fhicl::Atom<std::string> crvCoincidenceModuleLabel{Name("CrvCoincidenceModuleLabel"), Comment("CrvCoincidenceModuleLabel")};
      fhicl::Atom<std::string> crvCoincidenceMCModuleLabel{Name("CrvCoincidenceMCModuleLabel"), Comment("CrvCoincidenceMCModuleLabel")};
      fhicl::Atom<std::string> crvRecoPulseLabel{Name("CrvRecoPulseLabel"), Comment("CrvRecoPulseLabel")};
      fhicl::Atom<std::string> crvWaveformsModuleLabel{ Name("CrvWaveformsModuleLabel"), Comment("CrvWaveformsModuleLabel")};
      fhicl::Atom<std::string> crvDigiModuleLabel{ Name("CrvDigiModuleLabel"), Comment("CrvDigiModuleLabel")};
      fhicl::Table<SimParticleTimeOffset::Config> timeOffsets{ Name("TimeOffsets"), Comment("Time maps") };
      fhicl::Atom<art::InputTag> meanPBItag{Name("BeamIntensity"), Comment("Tag for BeamIntensity"), art::InputTag()};
      fhicl::Atom<double> deadTimeWindow{Name("deadTimeWindowEndMargin"),125};
      fhicl::Atom<double> signalWindowStartTime{Name("SignalWindowStartTime"),700};
      fhicl::Atom<double> signalWindowEndTime{Name("SignalWindowEndTime"),1705};
      //      fhicl::Atom<bool> ignoreUE{Name("IgnoreUE"),true};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit CRVAnalysisTree(const Parameters& conf);
    virtual ~CRVAnalysisTree() { }
    void beginJob();
    void beginSubRun(const art::SubRun & subrun ) override;
    void endJob();
    //    void beginSubRun(const art::SubRun & subrun ) override;
    void analyze(const art::Event& e);

  private:
    Config _conf;
    TTree* _crvana;
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

    EventInfo _einfo;
    double _deadtime1;
    double _deadtime2;
    double _fakeT0;
    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandExponential _randExpo;
    double meanMuonLife_;
    CLHEP::RandFlat _randFlat;

    void fillEventInfo(const art::Event& event);
    double getDeadTime(std::vector<CrvHitInfoReco> _crvinfo, bool ignoreUE);
    void resetBranches();
    SimParticleTimeOffset _toff;

  };

  CRVAnalysisTree::CRVAnalysisTree(const Parameters& conf):
    art::EDAnalyzer(conf),
    _conf(conf()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _randExpo(_engine),
    _randFlat(_engine),
    _toff(conf().timeOffsets())
  {

    GlobalConstantsHandle<PhysicsParams> phyPar;
    meanMuonLife_ = phyPar->getDecayTime();
    std::cout<<"Muon lifetime initialized to = "<<meanMuonLife_<<std::endl;
  }

  void CRVAnalysisTree::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _crvana = tfs->make<TTree>("crvTree", "CRV Info Tree");
    _crvana->Branch("evtinfo.",&_einfo,EventInfo::leafnames().c_str());
    _crvana->Branch("fakeT0",&_fakeT0,"fakeT0/D");

    // CRV info
    _crvana->Branch("deadtime1",&_deadtime1,"deadtime1/D");
    _crvana->Branch("deadtime2",&_deadtime2,"deadtime2/D");
    _crvana->Branch("crvinfo",&_crvinfo);
    _crvana->Branch("crvsummary",&_crvsummary);
    _crvana->Branch("bestcrv",&_bestcrv,"bestcrv/I");
    _crvana->Branch("crvpulseinfo",&_crvpulseinfo);
    _crvana->Branch("crvwaveforminfo",&_crvwaveforminfo);
    if(_conf.fillmc()){
      _crvana->Branch("crvinfomc",&_crvinfomc);
      _crvana->Branch("crvsummarymc",&_crvsummarymc);
      _crvana->Branch("crvinfomcplane",&_crvinfomcplane);
      _crvana->Branch("crvpulseinfomc",&_crvpulseinfomc);
    }

  }





  void CRVAnalysisTree::endJob(){}
  void CRVAnalysisTree::beginSubRun(const art::SubRun & subrun ) {}

  void CRVAnalysisTree::analyze(const art::Event& event){

    if(_conf.fillmc())
      _toff.updateMap(event);
    resetBranches();
    // fill basic event information
    _einfo._eventid = event.event();
    _einfo._runid = event.run();
    _einfo._subrunid = event.subRun();

    // mean number of protons on target
    _einfo._nprotons = -999;

    CRVAnalysis::FillCrvHitInfoCollections(_conf.crvCoincidenceModuleLabel(), _conf.crvCoincidenceMCModuleLabel(),
                                           _conf.crvRecoPulseLabel(), _conf.crvStepPointMCLabel(), _conf.simParticleLabel(), _conf.mcTrajectoryLabel(), event,
                                           _crvinfo, _crvinfomc, _crvsummary, _crvsummarymc, _crvinfomcplane, _conf.crvPlaneY());

    CRVAnalysis::FillCrvPulseInfoCollections(_conf.crvRecoPulseLabel(), _conf.crvWaveformsModuleLabel(), _conf.crvDigiModuleLabel(),
                                             _toff, event, _crvpulseinfo, _crvpulseinfomc, _crvwaveforminfo);
   
    //Calculate deadtime
    _deadtime1 = getDeadTime(_crvinfo, false);
    _deadtime2 = getDeadTime(_crvinfo, true);
    _crvana->Fill();
  }

  void CRVAnalysisTree::resetBranches() {
    // clear vectors
    _crvinfo.clear();
    _crvinfomc.clear();
    _crvinfomcplane.clear();
    _crvpulseinfo.clear();
    _crvwaveforminfo.clear();
    _crvpulseinfomc.clear();

    _einfo.reset();
    _deadtime1 = 0;
    _deadtime2 = 0;
  }

  double CRVAnalysisTree::getDeadTime(std::vector<CrvHitInfoReco> _crvinfo, bool ignoreUE) {
    double deadcounter = 0;
    const int nsteps = 10000;
    const double step_size = (_conf.signalWindowEndTime() - _conf.signalWindowStartTime())/nsteps;

    for (int t=0; t<nsteps; t++){

      double fakeT0 = 9999;
      while (fakeT0>_conf.signalWindowEndTime())
	fakeT0 = _conf.signalWindowStartTime() + _randExpo.fire(meanMuonLife_);
      _fakeT0 = fakeT0;

      bool vetoed = false;
      for(auto i : _crvinfo){
	if(ignoreUE && (i._crvSectorType == 4 || i._crvSectorType == 5))
	  continue;
	double endVeto = i._timeWindowEnd + _conf.deadTimeWindow();
	double startVeto = i._timeWindowStart;
	if( fakeT0 > startVeto && fakeT0 < endVeto ){
	  vetoed = true;
	  break;
	}
      }
      if(vetoed)
	deadcounter += step_size;
    }
    return deadcounter/(step_size*nsteps);
  }


}


using mu2e::CRVAnalysisTree;
DEFINE_ART_MODULE(CRVAnalysisTree);
