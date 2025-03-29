//
// An EDAnalyzer module that reads the Trigger Info
//
// Original author G. Pezzullo
//

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
// #include "canvas/Utilities/InputTag.h"

#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/LoopHelix.hh"

//Services
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

//Data products
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

//Gen-level info
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

//Utilities
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"

//ROOT
#include "TH1.h"
#include "TH2.h"

#include <cmath>
#include <string>
#include <sstream>
#include <vector>

namespace mu2e {

  class ReadTriggerInfo : public art::EDAnalyzer {

  public:
    enum {
      kNTrigInfo     = 1000,
      kNTrackTrig    = 1000,
      kNTrackTrigVar = 50,
      kNHelixTrig    = 1000,
      kNHelixTrigVar = 130,
      kNCaloCalib    = 1000,
      kNCaloCalibVar = 30,
      kNCaloOnly     = 1000,
      kNCaloOnlyVar  = 30,
      kNOcc          = 3001, //should be double the maximum N(track) due to 1 track + 1 helix, plus N(calo-only) plus 1 for general
      kNOccVar       = 100
    };

    // basic trigger information
    struct trigInfo_ {
      int           counts;
      int           exclusive_counts;
      std::string   label;

      trigInfo_ ():counts(0), exclusive_counts(0), label(""){}
    };

    // summary information about all triggers
    struct summaryInfoHist_ {
      TH1 *_hTrigInfo     [kNTrigInfo]; //1D trigger information
      TH1 *_hTrigModules  [kNTrigInfo]; //cut-flow on the module chain for each trigger
      TH2 *_h2DTrigInfo   [kNTrigInfo]; //2D trigger information
      TH1 *_hTrigBDW      [kNTrigInfo]; //trigger rates
      TH1 *_hTrigBits; //trigger path indices that are firing

      summaryInfoHist_() {
        // initialize each histogram to null
        _hTrigBits = nullptr;
        for (int i=0; i<kNTrigInfo; ++i) {
          _hTrigInfo   [i] = nullptr;
          _hTrigModules[i] = nullptr;
          _h2DTrigInfo [i] = nullptr;
          _hTrigBDW    [i] = nullptr;
        }
      }
    };

    // information about tracks corresponding to track-related triggers
    struct trackInfoHist_ {
      TH1 *_hTrkInfo [kNTrackTrig][kNTrackTrigVar]; //track histograms, one set of histograms per track-related trigger path

      trackInfoHist_() {
        for (int i=0; i<kNTrackTrig; ++i) {
          for (int j=0; j<kNTrackTrigVar; ++j) {
            _hTrkInfo  [i][j] = nullptr;
          }
        }
      }
    };

    // information about helices corresponding to helix-related triggers
    struct helixInfoHist_ {
      TH1 *_hHelInfo [kNHelixTrig][kNHelixTrigVar]; //helix histograms, one set of histograms per helix-related trigger path

      helixInfoHist_() {
        for (int i=0; i<kNHelixTrig; ++i) {
          for (int j=0; j<kNHelixTrigVar; ++j) {
            _hHelInfo  [i][j] = nullptr;
          }
        }
      }
    };

    // information about clusters corresponding to cluster-related triggers
    struct caloTrigSeedHist_ {
      TH1 *_hCaloOnlyInfo [kNCaloOnly][kNCaloOnlyVar]; //cluster histograms, one set of histograms per cluster-related trigger path

      caloTrigSeedHist_() {
        for (int i=0; i<kNCaloOnly; ++i) {
          for (int j=0; j<kNCaloOnlyVar; ++j) {
            _hCaloOnlyInfo  [i][j] = nullptr;
          }
        }
      }
    };

    // calo information for calo calibration triggers
    struct caloCalibrationHist_ {
      TH1 *_hCaloCalibInfo[kNCaloCalib][kNCaloCalibVar]; //calo histograms, one set of histograms per calorimeter calibration trigger path

      caloCalibrationHist_ () {
        for (int i=0; i<kNCaloCalib; ++i) {
          for (int j=0; j<kNCaloCalibVar; ++j) {
            _hCaloCalibInfo  [i][j] = nullptr;
          }
        }
      }
    };

    // occupancy information
    struct occupancyHist_ {
      TH1 *_hOccInfo  [kNOcc][kNOccVar];
      TH2 *_h2DOccInfo[kNOcc][kNOccVar];

      occupancyHist_ () {
        for (int i=0; i<kNOcc; ++i) {
          for (int j=0; j<kNOccVar; ++j) {
            _hOccInfo    [i][j] = nullptr;
            _h2DOccInfo  [i][j] = nullptr;
          }
        }
      }
    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               diagLevel          {Name("diagLevel"          ), Comment("turn tool on or off"                       ), 0 };
      fhicl::Atom<size_t>            nMaxTrig           {Name("nPathIDs"           ), Comment("List of the path Ids"                      ), 100};
      fhicl::Atom<size_t>            nTrackTrig         {Name("nTrackTriggers"     ), Comment("nTrackTriggers"                            )  };
      fhicl::Atom<size_t>            nCaloTrig          {Name("nCaloTriggers"      ), Comment("Number of calorimeter triggers"            )  };
      fhicl::Atom<size_t>            nCaloCalibTrig     {Name("nCaloCalibTriggers" ), Comment("Number of calorimeter calibration triggers")  };
      fhicl::Atom<art::InputTag>     sdTag              {Name("strawDigiCollection"), Comment("makeSD"                                    ), "makeSD"  };
      fhicl::Atom<art::InputTag>     chTag              {Name("comboHitCollection" ), Comment("TTmakeSH"                                  ), "TTmakeSH"};
      fhicl::Atom<art::InputTag>     cdTag              {Name("caloDigiCollection" ), Comment("CaloDigiFromShower"                        ), "CaloDigiFromShower"  };
      fhicl::Atom<art::InputTag>     genCountTag        {Name("genCount"           ), Comment("GenEventCount label"                       ), "genCounter" };
      fhicl::Atom<art::InputTag>     PBITag             {Name("PBITag"             ), Comment("ProtonBunchIntensity label"                ), "PBISim" };
      fhicl::Atom<float>             duty_cycle         {Name("dutyCycle"          ), Comment("Duty cycle"                                ), 0.3};
      fhicl::Atom<string>            processName        {Name("processName"        ), Comment("globalTrigger"                             ), ""  };
    };

    explicit ReadTriggerInfo(const art::EDAnalyzer::Table<Config>& config);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(const art::Run & run);
    virtual void beginSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);

    void     bookHistograms           ();
    void     bookTrigInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, summaryInfoHist_       &Hist);
    void     bookTrackInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, trackInfoHist_         &Hist);
    void     bookHelixInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, helixInfoHist_         &Hist);
    void     bookCaloTrigSeedInfoHist (art::ServiceHandle<art::TFileService> & Tfs, caloTrigSeedHist_      &Hist);
    void     bookCaloCalibInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist);
    void     bookOccupancyInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, occupancyHist_         &Hist);

    void     fillTrackTrigInfo        (int TrkTrigIndex  , const KalSeed*   KSeed, trackInfoHist_         &Hist);
    void     fillHelixTrigInfo        (int HelTrigIndex  , const HelixSeed* HSeed, helixInfoHist_         &Hist);
    void     fillCaloTrigInfo         (int Index         , const CaloCluster* HCl, caloTrigSeedHist_      &Hist);
    void     fillCaloCalibTrigInfo    (int Index         , const CaloCluster* HCl, caloCalibrationHist_   &Hist);
    void     fillOccupancyInfo        (int Index         , const StrawDigiCollection*SDCol, const CaloDigiCollection*CDCol, occupancyHist_   &Hist);

    bool     isTrackFilter            (const string& module);
    bool     isHelixFilter            (const string& module);
    bool     isCaloFilter             (const string& module);
    bool     isGlobalFilter           (const string& module);
    void     trigPathVal              (const int Index, TriggerResultsNavigator& trigNavig, summaryInfoHist_& Hist);
    void     findCorrelatedEvents     (std::vector<string>& VecLabels, double &NCorrelated);
    void     evalTriggerRate          ();

  private:
    int                       _diagLevel;
    size_t                    _nMaxTrig;
    int                       _nTrackTrig;
    int                       _nCaloTrig;
    int                       _nCaloCalibTrig;
    std::vector<std::string>  _trigPaths;
    art::InputTag             _trigAlgTag;
    art::InputTag             _sdTag;
    art::InputTag             _chTag;
    art::InputTag             _cdTag;
    art::InputTag             _genCountTag;
    art::InputTag             _PBITag;

    double                    _duty_cycle;
    string                    _processName;
    int                       _nProcess;
    double                    _bz0;
    double                    _nPOT;

    std::vector<trigInfo_>    _trigAll;
    std::vector<trigInfo_>    _trigFinal;
    std::vector<trigInfo_>    _trigCaloOnly;
    std::vector<trigInfo_>    _trigCaloCalib;
    std::vector<trigInfo_>    _trigTrack;
    std::vector<trigInfo_>    _trigHelix;
    std::vector<trigInfo_>    _trigEvtPS;

    summaryInfoHist_          _sumHist;
    trackInfoHist_            _trkHist;
    helixInfoHist_            _helHist;
    caloTrigSeedHist_         _caloTSeedHist;
    caloCalibrationHist_      _caloCalibHist;
    occupancyHist_            _occupancyHist;

    const mu2e::Tracker*      _tracker;
    const mu2e::ComboHitCollection*    _chcol;
    const art::Event*                  _event;
    const mu2e::HelixSeedCollection*   _hsCprCol;
    const mu2e::HelixSeedCollection*   _hsTprCol;

    float  _minPOT, _maxPOT;
    const int _genOccIndex;

    bool _useNGen; //use gen event count for normalization if available
  };

  ReadTriggerInfo::ReadTriggerInfo(const art::EDAnalyzer::Table<Config>& config):
    art::EDAnalyzer{config},
    _diagLevel           (config() .diagLevel()      ),
    _nMaxTrig            (config() .nMaxTrig()       ),
    _nTrackTrig          (config() .nTrackTrig()     ),
    _nCaloTrig           (config() .nCaloTrig()      ),
    _nCaloCalibTrig      (config() .nCaloCalibTrig() ),
    _sdTag               (config() .sdTag()          ),
    _chTag               (config() .chTag()          ),
    _cdTag               (config() .cdTag()          ),
    _genCountTag         (config() .genCountTag()    ),
    _PBITag              (config() .PBITag()         ),
    _duty_cycle          (config() .duty_cycle()     ),
    _processName         (config() .processName()    ),
    _nProcess            (0                          ),
    _minPOT              (0.                         ),
    _maxPOT              (3.5e8                      ),
    _genOccIndex         (2*_nTrackTrig + _nCaloTrig )
  {

    _trigAll.        resize(_nMaxTrig);
    _trigFinal.      resize(_nMaxTrig);
    _trigCaloOnly.   resize(_nMaxTrig);
    _trigCaloCalib.  resize(_nMaxTrig);
    _trigTrack.      resize(_nMaxTrig);
    _trigHelix.      resize(_nMaxTrig);
    _trigEvtPS.      resize(_nMaxTrig);

  }

  //--------------------------------------------------------------------------------//
  // Initialize all output histograms
  void ReadTriggerInfo::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    bookTrigInfoHist(tfs, _sumHist);
    bookTrackInfoHist(tfs, _trkHist);
    bookHelixInfoHist(tfs, _helHist);
    bookCaloTrigSeedInfoHist(tfs, _caloTSeedHist);
    bookCaloCalibInfoHist(tfs, _caloCalibHist);
    bookOccupancyInfoHist(tfs, _occupancyHist);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookTrigInfoHist(art::ServiceHandle<art::TFileService>& Tfs, summaryInfoHist_& Hist) {
    art::TFileDirectory trigInfoDir = Tfs->mkdir("trigInfo");
    Hist._hTrigBits      = trigInfoDir.make<TH1F>("hTrigInfo_bits"      , "Trigger bits"                               , 51, -0.5, 50.5);

    Hist._hTrigInfo[0]   = trigInfoDir.make<TH1F>("hTrigInfo_global"    , "Global Trigger rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[1]   = trigInfoDir.make<TH1F>("hTrigInfo_track"     , "Track Triggers rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[2]   = trigInfoDir.make<TH1F>("hTrigInfo_calo"      , "Calo-only Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[3]   = trigInfoDir.make<TH1F>("hTrigInfo_evtPS"     , "Event prescaler Trigger bits distribution"  , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[4]   = trigInfoDir.make<TH1F>("hTrigInfo_helix"     , "HelixSeed Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[5]   = trigInfoDir.make<TH1F>("hTrigInfo_caloCalib" , "Calo Calibration rejection"                 , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[6]   = trigInfoDir.make<TH1F>("hTrigInfo_final"     , "Global Trigger rejection of the paths"      , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

    Hist._hTrigInfo[10]  = trigInfoDir.make<TH1F>("hTrigInfo_unique_all", "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[11]  = trigInfoDir.make<TH1F>("hTrigInfo_unique"    , "Events found only by each Trig path"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

    Hist._hTrigInfo[15]  = trigInfoDir.make<TH1F>("hTrigInfo_paths"     , "Rejection of all the Trigger paths"         , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

    Hist._hTrigInfo[16]  = trigInfoDir.make<TH1F>("hTrksVsPOT","nOfflineTracks vs inst lum; p/pulse; nOfflieTracks",  1000, _minPOT, _maxPOT);
    Hist._hTrigInfo[17]  = trigInfoDir.make<TH1F>("hTrksTrigVsPOT","nOfflineTracks triggered vs inst lum; p/pulse; nOfflieTracks triggered",  1000, _minPOT, _maxPOT);
    Hist._hTrigInfo[18]  = trigInfoDir.make<TH1F>("hNormVsPOT","events vs inst lum; p/pulse; Entries",  1000, _minPOT, _maxPOT);
    Hist._hTrigInfo[19]  = trigInfoDir.make<TH1F>("hNPOT","nPOT; p/pulse; Entries",  1000, _minPOT, _maxPOT);
    Hist._hTrigInfo[20]  = trigInfoDir.make<TH1F>("hCprNorm","cpr Good Offline tracks;p/pulse; Entries",  1000, _minPOT, _maxPOT);
    Hist._hTrigInfo[21]  = trigInfoDir.make<TH1F>("hTprNorm","tpr Good Offline tracks;p/pulse; Entries",  1000, _minPOT, _maxPOT);

    // Trigger efficiency histograms
    Hist._hTrigInfo[30]  = trigInfoDir.make<TH1F>("hTrigInfo_eff"      , "Trigger efficiencies"                        , 100, 0, 100);

    Hist._h2DTrigInfo[0] = trigInfoDir.make<TH2F>("h2DTrigInfo_map_all" , "Trigger correlation map from all filters"   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._h2DTrigInfo[1] = trigInfoDir.make<TH2F>("h2DTrigInfo_map"     , "Trigger correlation map"                    , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

    //add a for loop for creating histograms showing the trigger path module filter path
    for(unsigned int i=0; i < _nMaxTrig; ++i) {
      Hist._hTrigModules[i] = trigInfoDir.make<TH1F>(Form("hTrigInfo_bk_%d",i)    , "breakdown histogram" , 30, 0, 30);
    }

    art::TFileDirectory trigBDWDir = Tfs->mkdir("trigBDW");

    Hist._hTrigBDW[0] = trigBDWDir.make<TH1F>("hTrigBDW_global"    , "Trigger bandwidth; ; rate [Hz]"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigBDW[1] = trigBDWDir.make<TH1F>("hTrigBDW_cumulative", "Cumulative Trigger bandwidth; ; rate [Hz]"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

  }
  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookTrackInfoHist(art::ServiceHandle<art::TFileService>& Tfs, trackInfoHist_& Hist) {
    for (int i=0; i<_nTrackTrig; ++i) {
      art::TFileDirectory trkInfoDir  = Tfs->mkdir(Form("trk_%i", i));
      Hist._hTrkInfo[i][0] = trkInfoDir.make<TH1F>(Form("hP_%i"     , i), "Track Momentum; p[MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][1] = trkInfoDir.make<TH1F>(Form("hPt_%i"    , i), "Track Pt; p_{t} [MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][2] = trkInfoDir.make<TH1F>(Form("hNSh_%i"   , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      Hist._hTrkInfo[i][3] = trkInfoDir.make<TH1F>(Form("hD0_%i"    , i), "Track impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hTrkInfo[i][4] = trkInfoDir.make<TH1F>(Form("hChi2d_%i" , i), "Track #chi^{2}/ndof;#chi^{2}/ndof", 100, 0, 50);
      Hist._hTrkInfo[i][5] = trkInfoDir.make<TH1F>(Form("hClE_%i"   , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
      Hist._hTrkInfo[i][6] = trkInfoDir.make<TH1F>(Form("hNLoops_%i", i), "Helix nLoops", 500, 0, 50);
      Hist._hTrkInfo[i][40] = trkInfoDir.make<TH1F>(Form("hNTrigTracks_%i" , i), "NTracks trigger matched; NTracks trigger matched", 11, -0.5, 10);
    }
  }
  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookHelixInfoHist(art::ServiceHandle<art::TFileService>& Tfs, helixInfoHist_& Hist) {
    for (int i=0; i<_nTrackTrig; ++i) {
      art::TFileDirectory helInfoDir  = Tfs->mkdir(Form("helix_%i", i));
      Hist._hHelInfo[i][0] = helInfoDir.make<TH1F>(Form("hP_%i"        , i), "Helix Momentum; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][1] = helInfoDir.make<TH1F>(Form("hPt_%i"       , i), "Helix Pt; p_{t} [MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][2] = helInfoDir.make<TH1F>(Form("hNSh_%i"      , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      Hist._hHelInfo[i][3] = helInfoDir.make<TH1F>(Form("hD0_%i"       , i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][4] = helInfoDir.make<TH1F>(Form("hChi2dXY_%i"  , i), "Helix #chi^{2}_{xy}/ndof;#chi^{2}_{xy}/ndof"      , 100, 0, 50);
      Hist._hHelInfo[i][5] = helInfoDir.make<TH1F>(Form("hChi2dZPhi_%i", i), "Helix #chi^{2}_{z#phi}/ndof;#chi^{2}_{z#phi}/ndof", 100, 0, 50);
      Hist._hHelInfo[i][6] = helInfoDir.make<TH1F>(Form("hClE_%i"      , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
      Hist._hHelInfo[i][7] = helInfoDir.make<TH1F>(Form("hLambda_%i"   , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);
      Hist._hHelInfo[i][8] = helInfoDir.make<TH1F>(Form("hNLoops_%i"   , i), "Helix nLoops; nLoops", 500, 0, 50);
      Hist._hHelInfo[i][9] = helInfoDir.make<TH1F>(Form("hHitRatio_%i" , i), "Helix hitRatio; NComboHits/nExpectedComboHits", 200, 0, 2);
      Hist._hHelInfo[i][120] = helInfoDir.make<TH1F>(Form("hNTrigHelixes_%i" , i), "NHelixes trigger matched; NHelixes trigger matched", 11, -0.5, 10);
    }
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookCaloTrigSeedInfoHist(art::ServiceHandle<art::TFileService>& Tfs, caloTrigSeedHist_& Hist) {
    for (int i=0; i<_nCaloTrig; ++i) {
      art::TFileDirectory caloInfoDir = Tfs->mkdir(Form("caloOnly_%i",i));
      Hist._hCaloOnlyInfo[i][0] = caloInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      Hist._hCaloOnlyInfo[i][1] = caloInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
    }
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookCaloCalibInfoHist(art::ServiceHandle<art::TFileService>& Tfs, caloCalibrationHist_& Hist) {
    for (int i=0; i<_nCaloCalibTrig; ++i) {
      art::TFileDirectory caloCalibInfoDir = Tfs->mkdir(Form("caloCalib_%i",i));
      Hist._hCaloCalibInfo[i][0] = caloCalibInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      Hist._hCaloCalibInfo[i][1] = caloCalibInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
    }
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::bookOccupancyInfoHist(art::ServiceHandle<art::TFileService>& Tfs, occupancyHist_& Hist) {

    for(int i=0; i<_nTrackTrig; ++i) {
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoTrk_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i) {
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoHel_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i) {
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoCaloTrig_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

    art::TFileDirectory occInfoDir = Tfs->mkdir("occInfoGeneral");
    Hist._hOccInfo  [_genOccIndex][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,_genOccIndex),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);
    Hist._h2DOccInfo[_genOccIndex][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,_genOccIndex),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    Hist._h2DOccInfo[_genOccIndex][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,_genOccIndex),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);

  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerInfo::isTrackFilter(const string& module) {
    return module.find("KSFilter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerInfo::isHelixFilter(const string& module) {
    return module.find("HSFilter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerInfo::isCaloFilter(const string& module) {
    return module.find("calo") != std::string::npos && module.find("Filter") != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  bool ReadTriggerInfo::isGlobalFilter(const string& module) {
    return (module.find("HSFilter")!= std::string::npos && module.find("CosmicHelix")!= std::string::npos) ||
      module.find("caloPhotonFilter")!= std::string::npos ||
      module.find("caloMVANNCEFilter")!= std::string::npos ||
      module.find("TSFilter")       != std::string::npos;
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::trigPathVal(const int Index, TriggerResultsNavigator& trigNavig, summaryInfoHist_& Hist) {
    const std::string path = trigNavig.getTrigPathName(Index);
    const unsigned lastModule = trigNavig.indexLastModule(path);
    auto h = Hist._hTrigModules[Index];
    if(!h) throw cet::exception("BADCONFIG") << __func__ << ": Trigger path index " << Index << " is out of bounds for initialized histograms\n";
    if(h->GetEntries() == 0) h->SetTitle(path.c_str());

    // cut-flow along the trigger path
    for (unsigned i = 0; i <= lastModule ; ++i) {
      const int bin = h->FindBin(i);
      //set the bin label to the module name if available
      if(h->GetBinContent(bin) <= 0. && trigNavig.triggerModules(path).size() > i) h->GetXaxis()->SetBinLabel(bin, trigNavig.triggerModules(path)[i].c_str());
      h->Fill(i);
    }
  }


  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::beginJob() {

    // Validate that the input collection sizes don't exceed the maxima
    if(_nMaxTrig > kNTrigInfo) throw cet::exception("BADCONFIG") << "Number of triggers assumed " << _nMaxTrig << " is greater than the maximum " << kNTrigInfo << "\n";
    if(_nTrackTrig > kNTrackTrig) throw cet::exception("BADCONFIG") << "Number of track triggers assumed " << _nTrackTrig << " is greater than the maximum " << kNTrackTrig << "\n";
    if(_nCaloTrig > kNCaloOnly) throw cet::exception("BADCONFIG") << "Number of calo-only triggers assumed " << _nCaloTrig << " is greater than the maximum " << kNCaloOnly << "\n";
    if(_nCaloCalibTrig > kNCaloCalib) throw cet::exception("BADCONFIG") << "Number of calo calibration triggers assumed " << _nCaloCalibTrig << " is greater than the maximum " << kNCaloCalib << "\n";

    // Initialize the output histograms
    bookHistograms();
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::endJob() {
    if(_nProcess <= 0) {
      if(_diagLevel > 0) printf("[ReadTriggerInfo::%s] Setting N(processed) from %i to 1\n", __func__, _nProcess);
      _nProcess = 1;
    }

    //////////////////////////////////////////
    // Set histograms' titles

    //Helix
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Setting helix histogram titles\n", __func__);
    for (int i=0; i<_nTrackTrig; ++i) {
      for (int j=0; j<kNHelixTrigVar; ++j) {
        if (_helHist._hHelInfo[i][j] == nullptr)    continue;
        std::string title = _trigHelix[i].label +": "+ _helHist._hHelInfo[i][j]->GetTitle();
        _helHist._hHelInfo[i][j]->SetTitle(title.c_str());
      }
    }

    //Tracks
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Setting track histogram titles\n", __func__);
    for (int i=0; i<_nTrackTrig; ++i) {
      for (int j=0; j<kNTrackTrigVar; ++j) {
        if (_trkHist._hTrkInfo[i][j] == nullptr)    continue;
        std::string title = _trigTrack[i].label +": "+ _trkHist._hTrkInfo[i][j]->GetTitle();
        _trkHist._hTrkInfo[i][j]->SetTitle(title.c_str());
      }
    }


    //occupancy
    //tracks
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Setting track occupancy histogram titles\n", __func__);
    for (int i=0; i<_nTrackTrig; ++i) {
      for (int j=0; j<kNOccVar; ++j) {
        if (_occupancyHist._hOccInfo[i][j] != nullptr)    {
          std::string title = _trigTrack[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != nullptr)    {
          std::string title = _trigTrack[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }
    //helix
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Setting helix occupancy histogram titles\n", __func__);
    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i) {
      for (int j=0; j<kNOccVar; ++j) {
        if (_occupancyHist._hOccInfo[i][j] != nullptr)    {
          std::string title = _trigHelix[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != nullptr)    {
          std::string title = _trigHelix[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }
    //calo trig
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Setting calo histogram titles\n", __func__);
    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i) {
      for (int j=0; j<kNOccVar; ++j) {
        if (_occupancyHist._hOccInfo[i][j] != nullptr)    {
          std::string title = _trigCaloOnly[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != nullptr)    {
          std::string title = _trigCaloOnly[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }

    //////////////////////////////////////////
    // Fill the histograms

    int indexTrigInfo11(0);
    if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] Filling histograms\n", __func__);
    for(int bin = 0; bin <= _sumHist._hTrigInfo[30]->GetNbinsX()+1; ++bin)
      _sumHist._hTrigInfo[30]->SetBinContent(bin, _sumHist._hTrigInfo[30]->GetBinContent(bin)/_nProcess);

    for (size_t i=0; i<_trigAll.size(); ++i ) {
      _sumHist._hTrigInfo  [0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      _sumHist._h2DTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());

      if (_trigAll[i].counts > 0) {
        _sumHist._hTrigInfo[0]->SetBinContent(i+1, _nProcess*1.f/_trigAll[i].counts);
        for (size_t j=0; j<_trigAll.size(); ++j ) {
          _sumHist._h2DTrigInfo[0]->GetYaxis()->SetBinLabel(j+1, _trigAll[j].label.c_str());
        }
      }

      _sumHist._hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigTrack[i].label.c_str());
      if (_trigTrack[i].counts > 0) _sumHist._hTrigInfo[1]->SetBinContent(i+1, _nProcess*1.f/_trigTrack[i].counts);

      _sumHist._hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigCaloOnly[i].label.c_str());
      if (_trigCaloOnly[i].counts > 0) _sumHist._hTrigInfo[2]->SetBinContent(i+1, _nProcess*1.f/_trigCaloOnly[i].counts);

      _sumHist._hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigEvtPS[i].label.c_str());
      if (_trigEvtPS[i].counts > 0) _sumHist._hTrigInfo[3]->SetBinContent(i+1, _trigEvtPS[i].counts);

      _sumHist._hTrigInfo[4]->GetXaxis()->SetBinLabel(i+1, _trigHelix[i].label.c_str());
      if (_trigHelix[i].counts > 0) _sumHist._hTrigInfo[4]->SetBinContent(i+1, _trigHelix[i].counts);

      _sumHist._hTrigInfo[5]->GetXaxis()->SetBinLabel(i+1, _trigCaloCalib[i].label.c_str());
      if (_trigCaloCalib[i].counts > 0) _sumHist._hTrigInfo[5]->SetBinContent(i+1, _trigCaloCalib[i].counts);

      if (_trigFinal[i].counts > 0) {
        _sumHist._hTrigInfo  [6]->GetXaxis()->SetBinLabel(i+1, _trigFinal[i].label.c_str());
        _sumHist._hTrigInfo  [6]->SetBinContent(i+1, _nProcess*1.f/_trigFinal[i].counts);
      }

      //fill  the histograms that shows how many events were found exclusively by each trigger path
      _sumHist._hTrigInfo  [10]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      double    content_trigInfo11 = _sumHist._hTrigInfo [10]->GetBinContent(i+1);
      if (content_trigInfo11>0) {
        _sumHist._hTrigInfo  [11]->GetXaxis()->SetBinLabel(indexTrigInfo11 +1, _trigAll[i].label.c_str());
        _sumHist._hTrigInfo  [11]->SetBinContent(indexTrigInfo11 +1, content_trigInfo11);
        ++indexTrigInfo11;
      }
    }

    //now let's filter the 2D correlation histogram with only those that actually triggered at least one event
    int                nbinsx = _sumHist._h2DTrigInfo[0]->GetNbinsX();
    int                nbinsy = _sumHist._h2DTrigInfo[0]->GetNbinsY();
    std::vector<int>   binsToSkip;

    for (int i=0; i<nbinsx; ++i) {
      bool used(false);

      for (int j=0; j<nbinsy; ++j) {
        if (_sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1) > 0) {
          used = true;
          break;
        }
      }
      if (!used) binsToSkip.push_back(i);
    }

    int   index_x(0);
    for (int i=0; i<nbinsx; ++i) {
      int    counts = std::count(binsToSkip.begin(), binsToSkip.end(), i);
      if (counts >= 1)       continue;
      //set the label
      _sumHist._h2DTrigInfo[1]->GetXaxis()->SetBinLabel(index_x+1, _trigAll[i].label.c_str());

      int    index_y(0);

      for (int j=0; j<nbinsy; ++j) {
        counts = std::count(binsToSkip.begin(), binsToSkip.end(), j);
        if (counts >= 1)       continue;
        double  content =  _sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1);
        _sumHist._h2DTrigInfo[1]->SetBinContent(index_x+1, index_y+1, content);

        //set the label
        if (index_x == 0) {
          _sumHist._h2DTrigInfo[1]->GetYaxis()->SetBinLabel(index_y+1, _trigAll[j].label.c_str());
        }

        ++index_y;
      }
      ++index_x;
    }

    // now evaluate the bandwidth
    // NOTE: "evalTriggerrate" re-orders the vectors _trigFinal
    evalTriggerRate();
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::evalTriggerRate() {
    //order the array with the filter used at the end of each path
    std::sort(_trigFinal.begin(), _trigFinal.end(), [](const auto a, const auto b) {return a.counts < b.counts; });

    double    mean_mb_rate   = 1.;///(mbtime/CLHEP::s)*_duty_cycle;

    bool      isFirst(true);
    int       index(0);

    std::vector<string>    labels_by_rate;

    for (size_t i=0; i< _trigFinal.size(); ++i) {
      double  nEvents = (double)_trigFinal[i].counts;
      if ( nEvents <= 1e-3)                 continue;

      labels_by_rate.push_back(_trigFinal[i].label);

      double  eff   = nEvents*1.f/_nProcess;
      double  rate  = mean_mb_rate*eff;
      _sumHist._hTrigBDW[0]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _sumHist._hTrigBDW[1]->GetXaxis()->SetBinLabel(index+1, _trigFinal[i].label.c_str());
      _sumHist._hTrigBDW[0]->SetBinContent(index+1, rate);

      if (isFirst) {
        _sumHist._hTrigBDW[1]->SetBinContent(index+1, rate);
        //        rate_ref  = rate;
        isFirst = false;
      }else{
        double    nCorrelated(0);
        findCorrelatedEvents(labels_by_rate, nCorrelated);

        rate = _sumHist._hTrigBDW[1]->GetBinContent(index) + (nEvents-nCorrelated)*1.f/_nProcess*mean_mb_rate;
        _sumHist._hTrigBDW[1]->SetBinContent(index+1, rate);
      }

      ++index;
    }

  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::findCorrelatedEvents(std::vector<string>& VecLabels, double& NCorrelated) {

    NCorrelated = 0;

    const char* label_ref = VecLabels.at(VecLabels.size()-1).c_str();
    if (VecLabels.size()<2) return;

    //    char* label(0);

    //    int        nLabels = VecLabels.size() -1;
    int        nbins   = _sumHist._h2DTrigInfo[1]->GetNbinsX();
    for(int i=0; i<nbins; ++i) {
      const char* label = _sumHist._h2DTrigInfo[1]->GetXaxis()->GetBinLabel(i+1);
      if (std::strcmp(label_ref, label) != 0)      continue;

      //      for (int k=0; k<nLabels; ++k) {
      //      label_ref = VecLabels.at(k).c_str();
      const char* label_2 = VecLabels.at(VecLabels.size()-2).c_str();
      for (int j=0; j<nbins; ++j) {
        //if (j == i)      break;
        label =   _sumHist._h2DTrigInfo[1]->GetYaxis()->GetBinLabel(j+1);
        if (std::strcmp(label_2, label) != 0)        continue;
        NCorrelated += _sumHist._h2DTrigInfo[1]->GetBinContent(i+1, j+1);
      }
      //}
      break;
    }

  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::beginRun(const art::Run & run) {
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker  = th.get();
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::beginSubRun(const art::SubRun& sr) {
    art::Handle<mu2e::GenEventCount> genCountH;
    sr.getByLabel(_genCountTag, genCountH);
    _useNGen = genCountH.isValid();
    if(_useNGen) _nProcess += genCountH.product()->count();
    if(_diagLevel > 0) printf("[ReadTriggerInfo::%s] use gen count product flag = %o\n", __func__, _useNGen);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::analyze(const art::Event& event) {

    if(!_useNGen) _nProcess += 1;

    //get the number of POT
    art::Handle<mu2e::ProtonBunchIntensity> PBIH;
    event.getByLabel(_PBITag, PBIH);
    if(PBIH.isValid()) _nPOT = PBIH.product()->intensity();
    else _nPOT  = -1.;
    if(_diagLevel > 1) printf("[ReadTriggerInfo::%s] N(POT) = %.2f\n", __func__, _nPOT);
    _sumHist._hTrigInfo[19]->Fill(_nPOT);

    //get the TriggerResult
    std::ostringstream oss;
    oss << "TriggerResults::"<<_processName;
    art::InputTag const tag{oss.str()};
    if(_diagLevel > 1) printf("[ReadTriggerInfo::%s] Trigger results tag = %s\n", __func__, oss.str().c_str());
    auto const trigResultsH   = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults*trigResults = trigResultsH.product();
    TriggerResultsNavigator   trigNavig(trigResults);

    //now set the value of _trigPaths
    _trigPaths = trigNavig.getTrigPaths();

    //initialize bin labels
    if(_sumHist._hTrigInfo[30]->Integral() <= 0.) {
      for (unsigned int i=0; i< trigNavig.getTrigPaths().size(); ++i) {
        const std::string path = trigNavig.getTrigPathName(i);
        _sumHist._hTrigInfo[30]->GetXaxis()->SetBinLabel(_sumHist._hTrigInfo[30]->FindBin(i), path.c_str());
      }
    }

    //fill the histogram with the accepted trigger bits
    for (unsigned int i=0; i< trigNavig.getTrigPaths().size(); ++i) {
      const std::string path = trigNavig.getTrigPathName(i);
      if(trigNavig.accepted(path)) {
        _sumHist._hTrigBits->Fill(trigNavig.findTrigPath(path));
        _sumHist._hTrigInfo[15]->Fill(trigNavig.findTrigPathID(path)); //accepted path IDs
        _sumHist._hTrigInfo[30]->Fill(i);
      }
    }

    //get the helix collections
    art::Handle<mu2e::HelixSeedCollection>  hsCprH;
    if (hsCprH.isValid()) {
      _hsCprCol = hsCprH.product();
    }else {
      _hsCprCol = nullptr;
    }

    art::Handle<mu2e::HelixSeedCollection>  hsTprH;
    if (hsTprH.isValid()) {
      _hsTprCol = hsTprH.product();
    }else {
      _hsTprCol = nullptr;
    }

    //get the StrawDigi Collection
    art::Handle<mu2e::StrawDigiCollection> sdH;
    event.getByLabel(_sdTag, sdH);
    const StrawDigiCollection* sdCol(0);
    if (sdH.isValid()) {
      sdCol = sdH.product();
    }

    //get the ComboHitCollection
    art::Handle<mu2e::ComboHitCollection> chH;
    event.getByLabel(_chTag, chH);
    if (chH.isValid()) {
      _chcol = chH.product();
    }else {
      _chcol = nullptr;
    }

    //get the CaloDigi Collection
    art::Handle<mu2e::CaloDigiCollection> cdH;
    event.getByLabel(_cdTag, cdH);
    const CaloDigiCollection* cdCol(0);
    if (cdH.isValid()) {
      cdCol = cdH.product();
    }

    //fill the general occupancy histogram
    fillOccupancyInfo   (_genOccIndex, sdCol, cdCol, _occupancyHist);

    std::vector<int>   trigFlagAll_index, trigFlag_index;

    art::Handle<TriggerInfo> hTrigInfoH;
    const mu2e::TriggerInfo* trigInfo(nullptr);

    // indices for each trigger type, used in occupancy histograms
    int track_trig_index = 0;
    int helix_trig_index = 0;
    int calo_trig_index  = 0;

    // Loop through all trigger paths, storing trigger info
    for (unsigned i=0; i < _trigPaths.size(); ++i) {
      std::string pathName = trigNavig.getTrigPathName(i);
      std::string& path = pathName;
      if(_diagLevel > 1) printf("[ReadTriggerInfo::%s] : Checking path %s\n", __func__, path.c_str());

      // store information about the trigger path acceptance chain
      trigPathVal(i, trigNavig, _sumHist);

      if (trigNavig.accepted(path)) {
        if(_diagLevel > 0) printf("[ReadTriggerInfo::%s] : >>> Path %s is accepted\n", __func__, path.c_str());
        std::vector<std::string>      moduleNames = trigNavig.triggerModules(path);

        if(_diagLevel>0) {
          printf("[ReadTriggerInfo::analyze] moduleNames size = %lu\n", moduleNames.size());
        }

        for (size_t j=0; j<moduleNames.size(); ++j) {
          std::string  moduleLabel = moduleNames[j];
          if(_diagLevel>0) {
            if (j==0) printf("[ReadTriggerInfo::analyze]      name      \n");
            printf("[ReadTriggerInfo::analyze] %10s\n", moduleLabel.c_str());
          }
          int          index_all(i);//0);
          int          index(i);//0);
          bool         passed(false);
          size_t       nTrigObj(0);
          //fill the Global Trigger bits info
          _trigAll[index_all].label  = moduleLabel;

          event.getByLabel(moduleLabel, hTrigInfoH);
          if (hTrigInfoH.isValid()) {
            trigInfo = hTrigInfoH.product();
            if(_diagLevel > 2) printf("[ReadTriggerInfo::%s] : Trig info product has been found\n", __func__);
          }

          if (isHelixFilter(moduleLabel)) {
            if(_diagLevel > 3) printf("[ReadTriggerInfo::%s] : Helix-filter module (%s) found in path %s\n", __func__, moduleLabel.c_str(), path.c_str());
            _trigHelix[index].label  = moduleLabel;
            ++_trigHelix[index].counts;
            passed = true;
            nTrigObj=0;
            if(!trigInfo) throw cet::exception("BADINPUTS") << "Trigger info product not found before needed!\n";
            for (auto const hseed: trigInfo->helixes()) {
              if(hseed) {
                ++nTrigObj;
                fillHelixTrigInfo(index, hseed.get(), _helHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo(_nTrackTrig+helix_trig_index, sdCol, cdCol, _occupancyHist);
                }
              }
              ++helix_trig_index;
            }//end loop over the helix-collection
            _helHist._hHelInfo[i][120]->Fill(nTrigObj);

          } else if (isTrackFilter(moduleLabel)) {
            if(_diagLevel > 3) printf("[ReadTriggerInfo::%s] : Track-filter module (%s) found in path %s\n", __func__, moduleLabel.c_str(), path.c_str());
            _trigTrack[index].label  = moduleLabel;
            ++_trigTrack[index].counts;
            passed = true;
            nTrigObj=0;
            if(!trigInfo) throw cet::exception("BADINPUTS") << "Trigger info product not found before needed!\n";
            for (auto const kseed: trigInfo->tracks()) {
              if(kseed) {
                ++nTrigObj;
                fillTrackTrigInfo(index, kseed.get(), _trkHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo(track_trig_index, sdCol, cdCol, _occupancyHist);
                }
              }
              ++track_trig_index;
            }//end loop over the kaseed-collection
            _trkHist._hTrkInfo[i][40]->Fill(nTrigObj);
            trigFlag_index.push_back(index_all);

          } else if ( moduleLabel.find("EventPrescale") != std::string::npos) {
            _trigEvtPS[index].label  = moduleLabel;
            ++_trigEvtPS[index].counts;
          } else if ( moduleLabel.find("caloCalibCosmic") != std::string::npos) {
            _trigCaloCalib[index].label  = moduleLabel;
            ++_trigCaloCalib[index].counts;
            passed = false;
            nTrigObj=0;
            if(!trigInfo) throw cet::exception("BADINPUTS") << "Trigger info product not found before needed!\n";
            for (auto const cluster : trigInfo->caloClusters()) {
              if(cluster) {
                ++nTrigObj;
                fillCaloCalibTrigInfo(index, cluster.get(), _caloCalibHist);
              }
            }//end loop over the cluster-collection
            trigFlag_index.push_back(index_all);
          } else if (isCaloFilter(moduleLabel)) {
            if(_diagLevel > 3) printf("[ReadTriggerInfo::%s] : Calo-filter module (%s) found in path %s\n", __func__, moduleLabel.c_str(), path.c_str());
            _trigCaloOnly[index].label  = moduleLabel;
            ++_trigCaloOnly[index].counts;
            passed = true;
            nTrigObj=0;
            if(!trigInfo) throw cet::exception("BADINPUTS") << "Trigger info product not found before needed!\n";
            for (auto const clseed: trigInfo->caloClusters()) {
              if(clseed) {
                ++nTrigObj;
                fillCaloTrigInfo(calo_trig_index, clseed.get(), _caloTSeedHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo   (_nTrackTrig*2+calo_trig_index, sdCol, cdCol, _occupancyHist);
                }
              }
              ++calo_trig_index;
            }//end loop
            //_caloTSeedHist._hCaloOnlyInfo[i][20]->Fill(nTrigObj);
            trigFlag_index.push_back(index_all);
          }

          if(isGlobalFilter(moduleLabel)) {
            _trigFinal[index].label    = moduleLabel;
            ++_trigFinal[index].counts;
            ++_trigAll[index_all].counts;
            trigFlagAll_index.push_back(index_all);
          }
        }//end loop over the modules in a given trigger path
      }
    }

    //now fill the correlation matrix
    for (size_t i=0; i<trigFlagAll_index.size(); ++i) {
      for (size_t j=0; j<trigFlagAll_index.size(); ++j) {
        _sumHist._h2DTrigInfo[0]->Fill(trigFlagAll_index.at(i), trigFlagAll_index.at(j));
      }
    }

    if (trigFlagAll_index.size() == 1) _sumHist._hTrigInfo[10]->Fill(trigFlagAll_index.at(0));
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::fillTrackTrigInfo(int TrkTrigIndex, const KalSeed* KSeed, trackInfoHist_& Hist) {
    if(!Hist._hTrkInfo[TrkTrigIndex][0]) throw cet::exception("BADCONFIG") << __func__ << ": Track histogram index " << TrkTrigIndex << " not initialized\n";
    GlobalConstantsHandle<ParticleDataList> pdt;

    int                nsh = (int)KSeed->hits().size();
    // KalSegment const& fseg = KSeed->segments().front();
    auto& fseg = KSeed->segments().front();

    double     ndof  = std::max(1.0,nsh - 5.0);
    double     p     = fseg.mom();
    double     chi2d = KSeed->chisquared()/ndof;
    double     pt    = p*std::cos(std::atan(fseg.helix().tanDip()));
    double     d0    = fseg.helix().d0();
    double     clE(-1.);

    // loop hits over z
    double z_min = std::numeric_limits<double>::max();
    double z_max = std::numeric_limits<double>::lowest();
    for (const auto& finter : KSeed->intersections()) {
      if (finter.position3().z() < z_min) {
        z_min = finter.position3().z();
      }
      if (finter.position3().z() > z_max) {
        z_max = finter.position3().z();
      }
    }

    double lambda = fseg.loopHelix().lam();
    double pitch  = std::abs(lambda)*6.28;
    double nLoops = (z_max - z_min)/pitch;

    if (KSeed->caloCluster()) clE = KSeed->caloCluster()->energyDep();

    Hist._hTrkInfo[TrkTrigIndex][0]->Fill(p);
    Hist._hTrkInfo[TrkTrigIndex][1]->Fill(pt);
    Hist._hTrkInfo[TrkTrigIndex][2]->Fill(nsh);
    Hist._hTrkInfo[TrkTrigIndex][3]->Fill(d0);
    Hist._hTrkInfo[TrkTrigIndex][4]->Fill(chi2d);
    Hist._hTrkInfo[TrkTrigIndex][5]->Fill(clE);
    Hist._hTrkInfo[TrkTrigIndex][6]->Fill(nLoops);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::fillHelixTrigInfo(int HelTrigIndex, const HelixSeed* HSeed, helixInfoHist_& Hist) {
    if(!Hist._hHelInfo[HelTrigIndex][0]) throw cet::exception("BADCONFIG") << __func__ << ": Helix histogram index " << HelTrigIndex << " not initialized\n";
    GlobalConstantsHandle<ParticleDataList> pdt;
    HelixTool helTool(HSeed, _tracker);

    int        nch       = (int)HSeed->hits().size();
    int        nsh(0);
    for (int i=0; i<nch; ++i) {
      nsh += HSeed->hits().at(i).nStrawHits();
    }
    float      mm2MeV    = (3./10.)*_bz0;

    double     p         = HSeed->helix().momentum()*mm2MeV;
    double     chi2dZPhi = HSeed->helix().chi2dZPhi();
    double     chi2dXY   = HSeed->helix().chi2dXY();
    double     pt        = HSeed->helix().radius()*mm2MeV;
    double     d0        = HSeed->helix().rcent() - HSeed->helix().radius();
    double     clE(-1.);
    double     lambda    = fabs(HSeed->helix().lambda());
    if (HSeed->caloCluster()) clE = HSeed->caloCluster()->energyDep();
    double     nLoops    = helTool.nLoops();

    Hist._hHelInfo[HelTrigIndex][0]->Fill(p);
    Hist._hHelInfo[HelTrigIndex][1]->Fill(pt);
    Hist._hHelInfo[HelTrigIndex][2]->Fill(nsh);
    Hist._hHelInfo[HelTrigIndex][3]->Fill(d0);
    Hist._hHelInfo[HelTrigIndex][4]->Fill(chi2dXY);
    Hist._hHelInfo[HelTrigIndex][5]->Fill(chi2dZPhi);
    Hist._hHelInfo[HelTrigIndex][6]->Fill(clE);
    Hist._hHelInfo[HelTrigIndex][7]->Fill(lambda);
    Hist._hHelInfo[HelTrigIndex][8]->Fill(nLoops);
    Hist._hHelInfo[HelTrigIndex][9]->Fill(helTool.hitRatio());
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::fillCaloTrigInfo(int CaloIndex, const CaloCluster* HCl, caloTrigSeedHist_& Hist) {
    if(!Hist._hCaloOnlyInfo[CaloIndex][0]) throw cet::exception("BADCONFIG") << __func__ << ": Calo-only histogram index " << CaloIndex << " not initialized\n";
    int        clsize    = HCl->size();
    double     energy    = HCl->energyDep();

    Hist._hCaloOnlyInfo[CaloIndex][0]->Fill(energy);
    Hist._hCaloOnlyInfo[CaloIndex][1]->Fill(clsize);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster* HCl, caloCalibrationHist_& Hist) {
    if(!Hist._hCaloCalibInfo[ClCalibIndex][0]) throw cet::exception("BADCONFIG") << __func__ << ": Calo-calib histogram index " << ClCalibIndex << " not initialized\n";
    int        clsize    = HCl->size();
    double     energy    = HCl->energyDep();

    Hist._hCaloCalibInfo[ClCalibIndex][0]->Fill(energy);
    Hist._hCaloCalibInfo[ClCalibIndex][1]->Fill(clsize);
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::fillOccupancyInfo(int Index, const StrawDigiCollection* SDCol, const CaloDigiCollection* CDCol, occupancyHist_& Hist) {
    if(!Hist._hOccInfo[Index][0]) throw cet::exception("BADCONFIG") << __func__ << ": Occupancy histogram index " << Index << " not initialized\n";
    if (_nPOT < 0)          return;
    int   nSD(-1), nCD(-1);
    if (SDCol) nSD = SDCol->size();
    if (CDCol) nCD = CDCol->size();

    Hist._hOccInfo  [Index][0]->Fill(_nPOT);

    Hist._h2DOccInfo[Index][0]->Fill(_nPOT, nSD);
    Hist._h2DOccInfo[Index][1]->Fill(_nPOT, nCD);
  }

}

DEFINE_ART_MODULE(mu2e::ReadTriggerInfo)
