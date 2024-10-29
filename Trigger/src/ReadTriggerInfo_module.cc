//
// An EDAnalyzer module that reads the Trigger Info
//
// Original author G. Pezzullo
//

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
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

//Dataproducts
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/DataProducts/inc/GenVector.hh"

//MC dataproducts
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

//Utilities
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"

//ROOT
#include "TH1F.h"
#include "TH2F.h"

#include <cmath>
// #include <iostream>
#include <string>
#include <sstream>
// #include <map>
#include <vector>

namespace mu2e {

  class ReadTriggerInfo : public art::EDAnalyzer {

  public:

    enum {
      kNTrigInfo     = 40,
      kNTrackTrig    = 40,
      kNTrackTrigVar = 50,
      kNHelixTrig    = 40,
      kNHelixTrigVar = 130,
      kNCaloCalib    = 5,
      kNCaloCalibVar = 30,
      kNCaloOnly     = 5,
      kNCaloOnlyVar  = 30,
      kNOcc          = 100,
      kNOccVar       = 100
    };

    struct MCInfo {
      double  pMC;
      double  pTMC;
      double  pZMC;
      double  dpMC;
      double  dpTMC;
      double  dpZMC;
      double  pdg;
      double  origZ;
      double  origR;
      double  pdgM;
      double  lambda;
      double  d0;
      double  p;
    };

    struct  trigInfo_ {
      int           counts;
      int           exclusive_counts;
      std::string   label;

      trigInfo_ ():counts(0), exclusive_counts(0), label(""){}
    };

    struct  summaryInfoHist_  {
      TH1F *_hTrigInfo  [kNTrigInfo];
      TH2F *_h2DTrigInfo[kNTrigInfo];
      TH1F *_hTrigBDW   [kNTrigInfo];

      TH1F *_hTrigBits;

      summaryInfoHist_() {
        _hTrigBits = NULL;
        for (int i=0; i<kNTrigInfo; ++i){
          _hTrigInfo  [i] = NULL;
          _h2DTrigInfo[i] = NULL;
          _hTrigBDW   [i] = NULL;
        }
      }
    };

    struct  trackInfoHist_    {
      TH1F *_hTrkInfo [kNTrackTrig][kNTrackTrigVar];

      trackInfoHist_ (){
        for (int i=0; i<kNTrackTrig; ++i){
          for (int j=0; j<kNTrackTrigVar; ++j){
            _hTrkInfo  [i][j] = NULL;
          }
        }
      }
    };

    struct  helixInfoHist_    {
      TH1F *_hHelInfo [kNHelixTrig][kNHelixTrigVar];

      helixInfoHist_(){
        for (int i=0; i<kNHelixTrig; ++i){
          for (int j=0; j<kNHelixTrigVar; ++j){
            _hHelInfo  [i][j] = NULL;
          }
        }
      }
    };

    struct  caloTrigSeedHist_ {
      TH1F *_hCaloOnlyInfo [kNCaloOnly][kNCaloOnlyVar];

      caloTrigSeedHist_(){
        for (int i=0; i<kNCaloOnly; ++i){
          for (int j=0; j<kNCaloOnlyVar; ++j){
            _hCaloOnlyInfo  [i][j] = NULL;
          }
        }
      }
    };

    struct  caloCalibrationHist_ {
      TH1F *_hCaloCalibInfo[kNCaloCalib][kNCaloCalibVar];

      caloCalibrationHist_ (){
        for (int i=0; i<kNCaloCalib; ++i){
          for (int j=0; j<kNCaloCalibVar; ++j){
            _hCaloCalibInfo  [i][j] = NULL;
          }
        }
      }
   };

    struct  occupancyHist_       {
      TH1F *_hOccInfo  [kNOcc][kNOccVar];
      TH2F *_h2DOccInfo[kNOcc][kNOccVar];

      occupancyHist_ (){
        for (int i=0; i<kNOcc; ++i){
          for (int j=0; j<kNOccVar; ++j){
            _hOccInfo    [i][j] = NULL;
            _h2DOccInfo  [i][j] = NULL;
          }
        }
      }
    };

    explicit ReadTriggerInfo(fhicl::ParameterSet const& pset);
    virtual ~ReadTriggerInfo() { }

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void     bookHistograms           ();
    void     bookTrigInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, summaryInfoHist_       &Hist);
    void     bookTrackInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, trackInfoHist_         &Hist);
    void     bookHelixInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, helixInfoHist_         &Hist);
    void     bookCaloTrigSeedInfoHist (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist);
    void     bookCaloCalibInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist);
    void     bookOccupancyInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, occupancyHist_         &Hist);

    void     findTrigIndex            (std::vector<trigInfo_> &Vec, std::string &ModuleLabel, int &Index);
    void     fillTrackTrigInfo        (int TrkTrigIndex  , const KalSeed*   KSeed, trackInfoHist_         &Hist);
    void     fillHelixTrigInfo        (int HelTrigIndex  , const HelixSeed* HSeed, helixInfoHist_         &Hist);
    void     fillHelixTrigInfoAdd     (int HelTrigIndex  , int MCMotherIndex, const HelixSeed* HSeed, helixInfoHist_         &Hist, MCInfo &TMPMCInfo);
    //    void     fillCaloTrigSeedInfo     (int ClCalibIndex  , const CaloCluster* HCl, caloCalibrationHist_   &Hist);
    void     fillCaloCalibTrigInfo    (int Index         , const CaloCluster* HCl, caloCalibrationHist_   &Hist);
    void     fillOccupancyInfo        (int Index         , const StrawDigiCollection*SDCol, const CaloDigiCollection*CDCol, occupancyHist_   &Hist);

    void     findCorrelatedEvents (std::vector<string>& VecLabels, double &NCorrelated);
    void     evalTriggerRate      ();

    void     fillTrackEfficiencyHist(const mu2e::KalSeedCollection*     KsCol,
                                     const mu2e::TrkQualCollection*     TrkQualCol,
                                     const art::TriggerResults*         TrigResults,
                                     const mu2e::StepPointMCCollection* Steps);
    bool     goodTrkTanDip(const mu2e::KalSeed*Ks);
  private:

    int                       _diagLevel;
    size_t                    _nMaxTrig;
    int                       _nTrackTrig;
    int                       _nCaloTrig;
    int                       _nCaloCalibTrig;
    std::vector<std::string>  _trigPaths;
    art::InputTag             _trigAlgTag;
    art::InputTag             _sdMCTag;
    art::InputTag             _sdTag;
    art::InputTag             _chTag;
    art::InputTag             _cdTag;
    art::InputTag             _evtWeightTag;
    art::InputTag             _hsCprTag;
    art::InputTag             _hsTprTag;
    art::InputTag             _ksTag;
    art::InputTag             _trkQualTag;
    art::InputTag             _vdTag;

    double                    _duty_cycle;
    string                    _processName;
    std::vector<size_t>       _effBits;
    float                     _trkMinTanDip;
    float                     _trkMaxTanDip;
    float                     _trkMaxD0;
    float                     _trkMinMVA;
    float                     _nProcess;
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

    //    caloTrigSeedHist_         _caloTSeedHist;
    caloCalibrationHist_      _caloTSeedHist;
    caloCalibrationHist_      _caloCalibHist;
    occupancyHist_            _occupancyHist;

    const mu2e::Tracker*      _tracker;

    //the following pointer is needed to navigate the MC truth info of the strawHits
    const mu2e::StrawDigiMCCollection* _mcdigis;
    const mu2e::ComboHitCollection*    _chcol;
    const art::Event*                  _event;
    const mu2e::HelixSeedCollection*   _hsCprCol;
    const mu2e::HelixSeedCollection*   _hsTprCol;

    float  _minPOT, _maxPOT;
  };

  ReadTriggerInfo::ReadTriggerInfo(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _diagLevel     (pset.get<int>   ("diagLevel", 0)),
    _nMaxTrig      (pset.get<size_t>("nPathIDs", 200)),
    _nTrackTrig    (pset.get<size_t>("nTrackTriggers", 4)),
    _nCaloTrig     (pset.get<size_t>("nCaloTriggers", 4)),
    _nCaloCalibTrig(pset.get<size_t>("nCaloCalibTriggers", 4)),
    _trigPaths     (pset.get<std::vector<std::string>>("triggerPathsList")),
    _sdMCTag       (pset.get<art::InputTag>("strawDigiMCCollection", "compressDigiMCs")),
    _sdTag         (pset.get<art::InputTag>("strawDigiCollection"  , "makeSD")),
    _chTag         (pset.get<art::InputTag>("comboHitCollection"   , "TTmakeSH")),
    _cdTag         (pset.get<art::InputTag>("caloDigiCollection"   , "CaloDigiFromShower")),
    _evtWeightTag  (pset.get<art::InputTag>("protonBunchIntensity" , "PBISim")),
    _hsCprTag      (pset.get<art::InputTag>("cprHelixSeedCollection", "CalHelixFinderDe:Positive")), // , "KFFDeMHPar")),
    _hsTprTag      (pset.get<art::InputTag>("tprHelixSeedCollection", "HelixFinderDe:Positive")), // , "KFFDeMHPar")),
    _ksTag         (pset.get<art::InputTag>("kalSeedCollection"  , "KFFDeMHPar")),
    _trkQualTag    (pset.get<art::InputTag>("trackQualCollection", "TrkQualDeMHPar")),
    _vdTag         (pset.get<art::InputTag>("vdStepPoints","NOTNOW")), // , "compressDigiMCs:virtualdetector")),
    _duty_cycle    (pset.get<float> ("dutyCycle", 1.)),
    _processName   (pset.get<string> ("processName", "globalTrigger")),
    _effBits       (pset.get<std::vector<size_t>>("effBits", std::vector<size_t>{2,8})),
    _trkMinTanDip  (pset.get<float> ("trkMinTanDip", 0.5)),
    _trkMaxTanDip  (pset.get<float> ("trkMaxTanDip", 1.)),
    _trkMaxD0      (pset.get<float> ("trkMaxD0", 100.)),
    _trkMinMVA     (pset.get<float> ("trkMinMVA", 0.8)),
    _nProcess      (pset.get<float> ("nEventsProcessed", 1.)),
    _minPOT(1e6), _maxPOT(4e8)
  {
    _trigAll.      resize(_nMaxTrig);
    _trigFinal.    resize(_nMaxTrig);
    _trigCaloOnly. resize(_nMaxTrig);
    _trigCaloCalib.resize(_nMaxTrig);
    _trigTrack.    resize(_nMaxTrig);
    _trigHelix.    resize(_nMaxTrig);
    _trigEvtPS.    resize(_nMaxTrig);

  }

  void     ReadTriggerInfo::bookHistograms           (){
    art::ServiceHandle<art::TFileService> tfs;

    bookTrigInfoHist(tfs, _sumHist);

    bookTrackInfoHist(tfs, _trkHist);

    bookHelixInfoHist(tfs, _helHist);

    bookCaloTrigSeedInfoHist(tfs, _caloTSeedHist);

    bookCaloCalibInfoHist(tfs, _caloCalibHist);

    bookOccupancyInfoHist(tfs, _occupancyHist);
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookTrigInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, summaryInfoHist_       &Hist){
    art::TFileDirectory trigInfoDir = Tfs->mkdir("trigInfo");
    Hist._hTrigBits      = trigInfoDir.make<TH1F>("hTrigInfo_bits"      , "Trigger bits"                               , 51, -0.5, 50.5);
    for (size_t i=0; i< _trigPaths.size(); ++i){
      Hist._hTrigBits->GetXaxis()->SetBinLabel(i+1, _trigPaths[i].c_str());
    }

    Hist._hTrigInfo[0]   = trigInfoDir.make<TH1F>("hTrigInfo_global"    , "Global Trigger rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[1]   = trigInfoDir.make<TH1F>("hTrigInfo_track"     , "Calo-only Triggers rejection"               , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigInfo[2]   = trigInfoDir.make<TH1F>("hTrigInfo_calo"      , "Track Triggers rejection"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
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

    Hist._h2DTrigInfo[0] = trigInfoDir.make<TH2F>("h2DTrigInfo_map_all" , "Trigger correlation map from all filters"   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._h2DTrigInfo[1] = trigInfoDir.make<TH2F>("h2DTrigInfo_map"     , "Trigger correlation map"                    , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5), (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

    art::TFileDirectory trigBDWDir = Tfs->mkdir("trigBDW");

    Hist._hTrigBDW[0]   = trigBDWDir.make<TH1F>("hTrigBDW_global"    , "Trigger bandwidth; ; rate [Hz]"                   , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));
    Hist._hTrigBDW[1]   = trigBDWDir.make<TH1F>("hTrigBDW_cumulative", "Cumulative Trigger bandwidth; ; rate [Hz]"        , (_nMaxTrig+2), -0.5, (_nMaxTrig+1.5));

  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookTrackInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, trackInfoHist_         &Hist){
    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory trkInfoDir  = Tfs->mkdir(Form("trk_%i", i));
      Hist._hTrkInfo[i][0] = trkInfoDir.make<TH1F>(Form("hP_%i"     , i), "Track Momentum; p[MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][1] = trkInfoDir.make<TH1F>(Form("hPt_%i"    , i), "Track Pt; p_{t} [MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][2] = trkInfoDir.make<TH1F>(Form("hNSh_%i"   , i), "N-StrawHits; nStrawHits", 101, -0.5, 100.5);
      Hist._hTrkInfo[i][3] = trkInfoDir.make<TH1F>(Form("hD0_%i"    , i), "Track impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hTrkInfo[i][4] = trkInfoDir.make<TH1F>(Form("hChi2d_%i" , i), "Track #chi^{2}/ndof;#chi^{2}/ndof", 100, 0, 50);
      Hist._hTrkInfo[i][5] = trkInfoDir.make<TH1F>(Form("hClE_%i"   , i), "calorimeter Cluster energy; E [MeV]", 240, 0, 120);
      Hist._hTrkInfo[i][6] = trkInfoDir.make<TH1F>(Form("hNLoops_%i", i), "Helix nLoops", 500, 0, 50);

      Hist._hTrkInfo[i][10] = trkInfoDir.make<TH1F>(Form("hPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hTrkInfo[i][11] = trkInfoDir.make<TH1F>(Form("hPtMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hTrkInfo[i][12] = trkInfoDir.make<TH1F>(Form("hPzMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hTrkInfo[i][13] = trkInfoDir.make<TH1F>(Form("hDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{trk} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hTrkInfo[i][14] = trkInfoDir.make<TH1F>(Form("hDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{trk} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hTrkInfo[i][15] = trkInfoDir.make<TH1F>(Form("hDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{trk} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hTrkInfo[i][16] = trkInfoDir.make<TH1F>(Form("hPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hTrkInfo[i][17] = trkInfoDir.make<TH1F>(Form("hGenZ_%i" , i), "z origin; z-origin [mm]"              , 300,   0,   15000);
      Hist._hTrkInfo[i][18] = trkInfoDir.make<TH1F>(Form("hGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hTrkInfo[i][19] = trkInfoDir.make<TH1F>(Form("hPDGM_%i" , i), "PDG Mother Id; PdgId-mother"          , 2253,   -30.5,   2222.5);
      Hist._hTrkInfo[i][20] = trkInfoDir.make<TH1F>(Form("hEMC_%i"  , i), "MC Energy; E_{MC} [MeV]"              , 400,   0,   200);

      Hist._hTrkInfo[i][40] = trkInfoDir.make<TH1F>(Form("hNTrigTracks_%i" , i), "NTracks trigger matched; NTracks trigger matched", 11, -0.5, 10);

    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookHelixInfoHist        (art::ServiceHandle<art::TFileService> & Tfs, helixInfoHist_         &Hist){
    for (int i=0; i<_nTrackTrig; ++i){
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


      Hist._hHelInfo[i][10] = helInfoDir.make<TH1F>(Form("hPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][11] = helInfoDir.make<TH1F>(Form("hPtMC_%i", i), "MC Track Pt @ tracker front; p_{t} [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][12] = helInfoDir.make<TH1F>(Form("hPzMC_%i", i), "MC Track Pt @ tracker front; p_{z} [MeV/c]" , 800, -200, 200);
      Hist._hHelInfo[i][13] = helInfoDir.make<TH1F>(Form("hDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][14] = helInfoDir.make<TH1F>(Form("hDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][15] = helInfoDir.make<TH1F>(Form("hDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][16] = helInfoDir.make<TH1F>(Form("hPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][17] = helInfoDir.make<TH1F>(Form("hGenZ_%i" , i), "z origin; z-origin [mm]"              , 300,   0,   15000);
      Hist._hHelInfo[i][18] = helInfoDir.make<TH1F>(Form("hGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][19] = helInfoDir.make<TH1F>(Form("hPDGM_%i" , i), "PDG Mother Id; PdgId-mother"          , 2253,   -30.5,   2222.5);
      //      Hist._hHelInfo[i][20] = helInfoDir.make<TH1F>(Form("hEMC_%i"  , i), "MC Energy; E_{MC} [MeV]"              , 400,   0,   200);

      Hist._hHelInfo[i][20] = helInfoDir.make<TH1F>(Form("hMuMinusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][21] = helInfoDir.make<TH1F>(Form("hMuMinusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][22] = helInfoDir.make<TH1F>(Form("hMuMinusD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][23] = helInfoDir.make<TH1F>(Form("hMuMinusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][24] = helInfoDir.make<TH1F>(Form("hMuMinusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][25] = helInfoDir.make<TH1F>(Form("hMuMinusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][26] = helInfoDir.make<TH1F>(Form("hMuMinusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][27] = helInfoDir.make<TH1F>(Form("hMuMinusGenZ_%i" , i), "origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][28] = helInfoDir.make<TH1F>(Form("hMuMinusGenR_%i" , i), "r origin; r-origin [mm]",500,   0,   5000);
      Hist._hHelInfo[i][29] = helInfoDir.make<TH1F>(Form("hMuMinusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][30] = helInfoDir.make<TH1F>(Form("hMuPlusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][31] = helInfoDir.make<TH1F>(Form("hMuPlusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][32] = helInfoDir.make<TH1F>(Form("hMuPlusD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][33] = helInfoDir.make<TH1F>(Form("hMuPlusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][34] = helInfoDir.make<TH1F>(Form("hMuPlusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][35] = helInfoDir.make<TH1F>(Form("hMuPlusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][36] = helInfoDir.make<TH1F>(Form("hMuPlusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][37] = helInfoDir.make<TH1F>(Form("hMuPlusGenZ_%i" , i), "origin; z-origin [mm];", 300,   0,   15000);
      Hist._hHelInfo[i][38] = helInfoDir.make<TH1F>(Form("hMuPlusGenR_%i" , i), "r origin;  r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][39] = helInfoDir.make<TH1F>(Form("hMuPlusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][40] = helInfoDir.make<TH1F>(Form("hIPAMuPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][41] = helInfoDir.make<TH1F>(Form("hIPAMuP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][42] = helInfoDir.make<TH1F>(Form("hIPAMuD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][43] = helInfoDir.make<TH1F>(Form("hIPAMuDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][44] = helInfoDir.make<TH1F>(Form("hIPAMuDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][45] = helInfoDir.make<TH1F>(Form("hIPAMuDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][46] = helInfoDir.make<TH1F>(Form("hIPAMuPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][47] = helInfoDir.make<TH1F>(Form("hIPAMuGenZ_%i" , i), "z origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][48] = helInfoDir.make<TH1F>(Form("hIPAMuGenR_%i" , i), "r origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][49] = helInfoDir.make<TH1F>(Form("hIPAMuLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][50] = helInfoDir.make<TH1F>(Form("hGammaPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][51] = helInfoDir.make<TH1F>(Form("hGammaP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][52] = helInfoDir.make<TH1F>(Form("hGammaD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][53] = helInfoDir.make<TH1F>(Form("hGammaDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][54] = helInfoDir.make<TH1F>(Form("hGammaDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][55] = helInfoDir.make<TH1F>(Form("hGammaDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][56] = helInfoDir.make<TH1F>(Form("hGammaPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][57] = helInfoDir.make<TH1F>(Form("hGammaGenZ_%i" , i), "origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][58] = helInfoDir.make<TH1F>(Form("hGammaGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][59] = helInfoDir.make<TH1F>(Form("hGammaLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][60] = helInfoDir.make<TH1F>(Form("hProtonPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][61] = helInfoDir.make<TH1F>(Form("hProtonP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][62] = helInfoDir.make<TH1F>(Form("hProtonD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][63] = helInfoDir.make<TH1F>(Form("hProtonDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][64] = helInfoDir.make<TH1F>(Form("hProtonDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][65] = helInfoDir.make<TH1F>(Form("hProtonDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][66] = helInfoDir.make<TH1F>(Form("hProtonPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][67] = helInfoDir.make<TH1F>(Form("hProtonGenZ_%i" , i), "origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][68] = helInfoDir.make<TH1F>(Form("hProtonGenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][69] = helInfoDir.make<TH1F>(Form("hProtonLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][70] = helInfoDir.make<TH1F>(Form("hN0PMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][71] = helInfoDir.make<TH1F>(Form("hN0P_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][72] = helInfoDir.make<TH1F>(Form("hN0D0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][73] = helInfoDir.make<TH1F>(Form("hN0DP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][74] = helInfoDir.make<TH1F>(Form("hN0DPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][75] = helInfoDir.make<TH1F>(Form("hN0DPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][76] = helInfoDir.make<TH1F>(Form("hN0PDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][77] = helInfoDir.make<TH1F>(Form("hN0GenZ_%i" , i), "origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][78] = helInfoDir.make<TH1F>(Form("hN0GenR_%i" , i), "radial position origin; r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][79] = helInfoDir.make<TH1F>(Form("hN0Lambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][80] = helInfoDir.make<TH1F>(Form("hPiMinusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][81] = helInfoDir.make<TH1F>(Form("hPiMinusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][82] = helInfoDir.make<TH1F>(Form("hPiMinusD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][83] = helInfoDir.make<TH1F>(Form("hPiMinusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][84] = helInfoDir.make<TH1F>(Form("hPiMinusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][85] = helInfoDir.make<TH1F>(Form("hPiMinusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][86] = helInfoDir.make<TH1F>(Form("hPiMinusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][87] = helInfoDir.make<TH1F>(Form("hPiMinusGenZ_%i" , i), "origin; z-origin [mm]", 300,   0,   15000);
      Hist._hHelInfo[i][88] = helInfoDir.make<TH1F>(Form("hPiMinusGenR_%i" , i), "r origin; r-origin [mm]",500,   0,   5000);
      Hist._hHelInfo[i][89] = helInfoDir.make<TH1F>(Form("hPiMinusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][90] = helInfoDir.make<TH1F>(Form("hPiPlusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][91] = helInfoDir.make<TH1F>(Form("hPiPlusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 200);
      Hist._hHelInfo[i][92] = helInfoDir.make<TH1F>(Form("hPiPlusD0_%i", i), "Helix impact parameter; d0 [mm]", 801, -400.5, 400.5);
      Hist._hHelInfo[i][93] = helInfoDir.make<TH1F>(Form("hPiPlusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][94] = helInfoDir.make<TH1F>(Form("hPiPlusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][95] = helInfoDir.make<TH1F>(Form("hPiPlusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][96] = helInfoDir.make<TH1F>(Form("hPiPlusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][97] = helInfoDir.make<TH1F>(Form("hPiPlusGenZ_%i" , i), "origin; z-origin [mm];", 300,   0,   15000);
      Hist._hHelInfo[i][98] = helInfoDir.make<TH1F>(Form("hPiPlusGenR_%i" , i), "r origin;  r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][99] = helInfoDir.make<TH1F>(Form("hPiPlusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 500);

      Hist._hHelInfo[i][100] = helInfoDir.make<TH1F>(Form("hEMinusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][101] = helInfoDir.make<TH1F>(Form("hEMinusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 250);
      Hist._hHelInfo[i][102] = helInfoDir.make<TH1F>(Form("hEMinusD0_%i", i), "Helix impact parameter; d0 [mm]", 800, -800, 800);
      Hist._hHelInfo[i][103] = helInfoDir.make<TH1F>(Form("hEMinusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][104] = helInfoDir.make<TH1F>(Form("hEMinusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][105] = helInfoDir.make<TH1F>(Form("hEMinusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][106] = helInfoDir.make<TH1F>(Form("hEMinusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][107] = helInfoDir.make<TH1F>(Form("hEMinusGenZ_%i" , i), "origin; z-origin [mm];", 300,   0,   15000);
      Hist._hHelInfo[i][108] = helInfoDir.make<TH1F>(Form("hEMinusGenR_%i" , i), "r origin;  r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][109] = helInfoDir.make<TH1F>(Form("hEMinusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 600);

      Hist._hHelInfo[i][110] = helInfoDir.make<TH1F>(Form("hEPlusPMC_%i" , i), "MC Track Momentum @ tracker front; p[MeV/c]", 400, 0, 200);
      Hist._hHelInfo[i][111] = helInfoDir.make<TH1F>(Form("hEPlusP_%i", i), "Track P; p [MeV/c]" , 400, 0, 250);
      Hist._hHelInfo[i][112] = helInfoDir.make<TH1F>(Form("hEPlusD0_%i", i), "Helix impact parameter; d0 [mm]", 800, -800, 800);
      Hist._hHelInfo[i][113] = helInfoDir.make<TH1F>(Form("hEPlusDP_%i"  , i), "#Delta p @ tracker front; #Delta p = p_{hel} - p_{MC} [MeV/c]"     , 800, -200, 200);
      Hist._hHelInfo[i][114] = helInfoDir.make<TH1F>(Form("hEPlusDPt_%i"  , i), "#Delta pT @ tracker front; #Delta pT = pT_{hel} - pT_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][115] = helInfoDir.make<TH1F>(Form("hEPlusDPz_%i"  , i), "#Delta pZ @ tracker front; #Delta pZ = pZ_{hel} - pZ_{MC} [MeV/c]", 800, -200, 200);
      Hist._hHelInfo[i][116] = helInfoDir.make<TH1F>(Form("hEPlusPDG_%i"  , i), "PDG Id; PdgId"                        , 2253,   -30.5,   2222.5);
      Hist._hHelInfo[i][117] = helInfoDir.make<TH1F>(Form("hEPlusGenZ_%i" , i), "origin; z-origin [mm];", 300,   0,   15000);
      Hist._hHelInfo[i][118] = helInfoDir.make<TH1F>(Form("hEPlusGenR_%i" , i), "r origin;  r-origin [mm]", 500,   0,   5000);
      Hist._hHelInfo[i][119] = helInfoDir.make<TH1F>(Form("hEPlusLambda_%i" , i), "Helix #lambda=dz/d#phi; |#lambda| [mm/rad]", 500, 0, 600);

      Hist._hHelInfo[i][120] = helInfoDir.make<TH1F>(Form("hNTrigHelixes_%i" , i), "NHelixes trigger matched; NHelixes trigger matched", 11, -0.5, 10);
    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookCaloTrigSeedInfoHist (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_      &Hist){
    for (int i=0; i<_nCaloTrig; ++i){
      art::TFileDirectory caloInfoDir = Tfs->mkdir(Form("caloOnly_%i",i));
      Hist._hCaloCalibInfo[i][0] = caloInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      Hist._hCaloCalibInfo[i][1] = caloInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
      // Hist._hCaloOnlyInfo[i][0] = caloInfoDir.make<TH1F>(Form("hEPeak_%i"   , i), "peak energy; E[MeV]"        , 400, 0, 200);
      // Hist._hCaloOnlyInfo[i][1] = caloInfoDir.make<TH1F>(Form("hR1Max1_%i"   , i), "ring1 max; ring1max [MeV]" , 400, 0, 200);
      // Hist._hCaloOnlyInfo[i][2] = caloInfoDir.make<TH1F>(Form("hR1Max2_%i"   , i), "ring1 max; ring1max2 [MeV]", 400, 0, 200);
    }
  }
  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookCaloCalibInfoHist    (art::ServiceHandle<art::TFileService> & Tfs, caloCalibrationHist_   &Hist){
    for (int i=0; i<_nCaloCalibTrig; ++i){
      art::TFileDirectory caloCalibInfoDir = Tfs->mkdir(Form("caloCalib_%i",i));
      Hist._hCaloCalibInfo[i][0] = caloCalibInfoDir.make<TH1F>(Form("hE_%i"   , i), "Cluster energy; E[MeV]", 800, 0, 800);
      Hist._hCaloCalibInfo[i][1] = caloCalibInfoDir.make<TH1F>(Form("hN_%i"   , i), "Cluster size; nCrystalHits", 101, -0.5, 100.5);
    }
  }

  //--------------------------------------------------------------------------------//
  void     ReadTriggerInfo::bookOccupancyInfoHist         (art::ServiceHandle<art::TFileService> & Tfs, occupancyHist_       &Hist){

    for (int i=0; i<_nTrackTrig; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoTrk_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoHel_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i){
      art::TFileDirectory occInfoDir = Tfs->mkdir(Form("occInfoCaloTrig_%i", i));
      Hist._hOccInfo  [i][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,i),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

      Hist._h2DOccInfo[i][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,i),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
      Hist._h2DOccInfo[i][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,i),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
    }

     int    index_last = _nTrackTrig+_nCaloTrig;
     art::TFileDirectory occInfoDir = Tfs->mkdir("occInfoGeneral");
     Hist._hOccInfo  [index_last][0]  = occInfoDir.make<TH1F>(Form("hInstLum_%i"  ,index_last),"distrbution of instantaneous lum; p/#mu-bunch"  ,  1000, _minPOT, _maxPOT);

     Hist._h2DOccInfo[index_last][0]  = occInfoDir.make<TH2F>(Form("hNSDVsLum_%i" ,index_last),"inst lum vs nStrawDigi; p/#mu-bunch; nStrawDigi",  1000, _minPOT, _maxPOT, 5000, 0., 20000.);
     Hist._h2DOccInfo[index_last][1]  = occInfoDir.make<TH2F>(Form("hNCDVsLum_%i" ,index_last),"inst lum vs nCaloDigi; p/#mu-bunch; nCaloDigi"  ,  1000, _minPOT, _maxPOT, 5000, 0., 20000.);



  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::beginJob(){

    bookHistograms();
  }

  //--------------------------------------------------------------------------------//
  void ReadTriggerInfo::endJob(){

    //set hitograms' titles
    //Helix
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNHelixTrigVar; ++j){
        if (_helHist._hHelInfo[i][j] == NULL)    continue;
        string title = _trigHelix[i].label +": "+ _helHist._hHelInfo[i][j]->GetTitle();
        _helHist._hHelInfo[i][j]->SetTitle(title.c_str());
      }
    }

    //Tracks
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNTrackTrigVar; ++j){
        if (_trkHist._hTrkInfo[i][j] == NULL)    continue;
        string title = _trigTrack[i].label +": "+ _trkHist._hTrkInfo[i][j]->GetTitle();
        _trkHist._hTrkInfo[i][j]->SetTitle(title.c_str());
      }
    }


    //occupancy
    //tracks
    for (int i=0; i<_nTrackTrig; ++i){
      for (int j=0; j<kNOccVar; ++j){
        if (_occupancyHist._hOccInfo[i][j] != NULL)    {
          string title = _trigTrack[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
          string title = _trigTrack[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }
    //helix
    for (int i=_nTrackTrig; i<_nTrackTrig*2; ++i){
      for (int j=0; j<kNOccVar; ++j){
        if (_occupancyHist._hOccInfo[i][j] != NULL)    {
          string title = _trigHelix[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
          string title = _trigHelix[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }
    //calo trig
    for (int i=_nTrackTrig*2; i<_nTrackTrig*2+_nCaloTrig; ++i){
      for (int j=0; j<kNOccVar; ++j){
        if (_occupancyHist._hOccInfo[i][j] != NULL)    {
          string title = _trigCaloOnly[i].label +": "+ _occupancyHist._hOccInfo[i][j]->GetTitle();
          _occupancyHist._hOccInfo[i][j]->SetTitle(title.c_str());
        }
        if (_occupancyHist._h2DOccInfo[i][j] != NULL)    {
          string title = _trigCaloOnly[i].label +": "+ _occupancyHist._h2DOccInfo[i][j]->GetTitle();
          _occupancyHist._h2DOccInfo[i][j]->SetTitle(title.c_str());
        }
      }
    }

    int    indexTrigInfo11(0);
    //fill the histograms
    for (size_t i=0; i<_trigAll.size(); ++i ){
      _sumHist._hTrigInfo  [0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      _sumHist._h2DTrigInfo[0]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());

      if (_trigAll[i].counts > 0) {
        _sumHist._hTrigInfo[0]->SetBinContent(i+1, _nProcess/_trigAll[i].counts);
        for (size_t j=0; j<_trigAll.size(); ++j ){
          _sumHist._h2DTrigInfo[0]->GetYaxis()->SetBinLabel(j+1, _trigAll[j].label.c_str());
        }
      }

      _sumHist._hTrigInfo[1]->GetXaxis()->SetBinLabel(i+1, _trigTrack[i].label.c_str());
      if (_trigTrack[i].counts > 0) _sumHist._hTrigInfo[1]->SetBinContent(i+1, _nProcess/_trigTrack[i].counts);

      _sumHist._hTrigInfo[2]->GetXaxis()->SetBinLabel(i+1, _trigCaloOnly[i].label.c_str());
      if (_trigCaloOnly[i].counts > 0) _sumHist._hTrigInfo[2]->SetBinContent(i+1, _nProcess/_trigCaloOnly[i].counts);

      _sumHist._hTrigInfo[3]->GetXaxis()->SetBinLabel(i+1, _trigEvtPS[i].label.c_str());
      if (_trigEvtPS[i].counts > 0) _sumHist._hTrigInfo[3]->SetBinContent(i+1, _trigEvtPS[i].counts);

      _sumHist._hTrigInfo[4]->GetXaxis()->SetBinLabel(i+1, _trigHelix[i].label.c_str());
      if (_trigHelix[i].counts > 0) _sumHist._hTrigInfo[4]->SetBinContent(i+1, _trigHelix[i].counts);

      _sumHist._hTrigInfo[5]->GetXaxis()->SetBinLabel(i+1, _trigCaloCalib[i].label.c_str());
      if (_trigCaloCalib[i].counts > 0) _sumHist._hTrigInfo[5]->SetBinContent(i+1, _trigCaloCalib[i].counts);

      if (_trigFinal[i].counts > 0) {
        _sumHist._hTrigInfo  [6]->GetXaxis()->SetBinLabel(i+1, _trigFinal[i].label.c_str());
        _sumHist._hTrigInfo  [6]->SetBinContent(i+1, _nProcess/_trigFinal[i].counts);
      }

      //fill  the histograms that shows how many events were found exclusively by each trigger path
      _sumHist._hTrigInfo  [10]->GetXaxis()->SetBinLabel(i+1, _trigAll[i].label.c_str());
      double    content_trigInfo11 = _sumHist._hTrigInfo [10]->GetBinContent(i+1);
      if (content_trigInfo11>0){
        _sumHist._hTrigInfo  [11]->GetXaxis()->SetBinLabel(indexTrigInfo11 +1, _trigAll[i].label.c_str());
        _sumHist._hTrigInfo  [11]->SetBinContent(indexTrigInfo11 +1, content_trigInfo11);
        ++indexTrigInfo11;
      }

    }

    //now let's filter the 2D correlation histogram with only those that actually triggered at least one event
    int                nbinsx = _sumHist._h2DTrigInfo[0]->GetNbinsX();
    int                nbinsy = _sumHist._h2DTrigInfo[0]->GetNbinsY();
    std::vector<int>   binsToSkip;

    for (int i=0; i<nbinsx; ++i){
      bool used(false);

      for (int j=0; j<nbinsy; ++j){
        if (_sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1) > 0) {
          used = true;
          break;
        }
      }
      if (!used) binsToSkip.push_back(i);
    }

    int   index_x(0);
    for (int i=0; i<nbinsx; ++i){
      int    counts = std::count(binsToSkip.begin(), binsToSkip.end(), i);
      if (counts >= 1)       continue;
      //set the label
      _sumHist._h2DTrigInfo[1]->GetXaxis()->SetBinLabel(index_x+1, _trigAll[i].label.c_str());

      int    index_y(0);

      for (int j=0; j<nbinsy; ++j){
        counts = std::count(binsToSkip.begin(), binsToSkip.end(), j);
        if (counts >= 1)       continue;
        double  content =  _sumHist._h2DTrigInfo[0]->GetBinContent(i+1, j+1);
        _sumHist._h2DTrigInfo[1]->SetBinContent(index_x+1, index_y+1, content);

        //set the label
        if (index_x == 0){
          _sumHist._h2DTrigInfo[1]->GetYaxis()->SetBinLabel(index_y+1, _trigAll[j].label.c_str());
        }

        ++index_y;
      }
      ++index_x;
    }

    // now evaluate the bandwidth
    // NOTE: "evalTriggerrate" re-order the vectors _trigFinal
    evalTriggerRate();
  }

  //--------------------------------------------------------------------------------
  void   ReadTriggerInfo::evalTriggerRate        (){
    //order the array with the filter used at the end of each path
    std::sort(_trigFinal.begin(), _trigFinal.end(), [](const auto a, const auto b) {return a.counts < b.counts; });

    double    mean_mb_rate   = 1.;///(mbtime/CLHEP::s)*_duty_cycle;

    bool      isFirst(true);
    int       index(0);

    std::vector<string>    labels_by_rate;

    for (size_t i=0; i< _trigFinal.size(); ++i){
      double  nEvents = (double)_trigFinal[i].counts;
      if ( nEvents <= 1e-3)                 continue;

      labels_by_rate.push_back(_trigFinal[i].label);

      double  eff   = nEvents/_nProcess;
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

        rate = _sumHist._hTrigBDW[1]->GetBinContent(index) + (nEvents-nCorrelated)/(double)_nProcess*mean_mb_rate;
              _sumHist._hTrigBDW[1]->SetBinContent(index+1, rate);
      }

      ++index;
    }

  }

  void   ReadTriggerInfo::findCorrelatedEvents(std::vector<string>& VecLabels, double &NCorrelated){

    NCorrelated = 0;

    const char* label_ref = VecLabels.at(VecLabels.size()-1).c_str();
    if (VecLabels.size()<2) return;

    //    char* label(0);

    //    int        nLabels = VecLabels.size() -1;
    int        nbins   = _sumHist._h2DTrigInfo[1]->GetNbinsX();
    for(int i=0; i<nbins; ++i){
      const char* label = _sumHist._h2DTrigInfo[1]->GetXaxis()->GetBinLabel(i+1);
      if (std::strcmp(label_ref, label) != 0)      continue;

      //      for (int k=0; k<nLabels; ++k){
      //      label_ref = VecLabels.at(k).c_str();
      const char* label_2 = VecLabels.at(VecLabels.size()-2).c_str();
      for (int j=0; j<nbins; ++j){
        //if (j == i)      break;
        label =   _sumHist._h2DTrigInfo[1]->GetYaxis()->GetBinLabel(j+1);
        if (std::strcmp(label_2, label) != 0)        continue;
        NCorrelated += _sumHist._h2DTrigInfo[1]->GetBinContent(i+1, j+1);
      }
      //}
      break;
    }

  }

  //================================================================
  void   ReadTriggerInfo::beginRun(const art::Run & run){
    // get bfield
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker  = th.get();
  }

  void ReadTriggerInfo::endSubRun(const art::SubRun& sr){}

  bool ReadTriggerInfo::goodTrkTanDip(const mu2e::KalSeed*Ks){
    const mu2e::KalSegment* kSeg = &(Ks->segments().at(0));
    float tanDip = kSeg->helix().tanDip();
    if ( (tanDip > _trkMinTanDip) && (tanDip< _trkMaxTanDip) ) return true;

    return false;
  }

  void ReadTriggerInfo::fillTrackEfficiencyHist(const mu2e::KalSeedCollection* KsCol,
                                                const mu2e::TrkQualCollection* TrkQualCol,
                                                const art::TriggerResults*     TrigResults,
                                                const mu2e::StepPointMCCollection* Steps){

    float pCEleCuts[2]={103., 105.};
    float pCPosCuts[2]={90.5, 92.5};

    const mu2e::HelixSeed* hs(0);

    if (KsCol != NULL){
      const mu2e::KalSeed*ks(0);
      for (size_t i=0; i<KsCol->size(); ++i){
        ks = &KsCol->at(i);
        const mu2e::KalSegment* kSeg = &(ks->segments().at(0));
        float p = kSeg->mom();

        if (ks->particle() == PDGCode::e_minus){
            if ((p<pCEleCuts[0]) || (p>pCEleCuts[1]) )
              continue;
        } else         if  (ks->particle() == PDGCode::e_plus){
          if ((p<pCPosCuts[0]) || (p>pCPosCuts[1]) )
            continue;
        }
        if ( (TrkQualCol->at(i).MVAOutput()<_trkMinMVA) ||
             (!goodTrkTanDip(ks)) || // [0.5, 1.0]
             (std::abs(kSeg->helix().d0()) > _trkMaxD0) ) continue;
        _sumHist._hTrigInfo[16]->Fill(_nPOT);

        float ks_t0 = ks->t0()._t0;
        float t0Toll(80.);
        bool  hasCprHelix(false), hasTprHelix(false);
        for (size_t i=0; i<_hsCprCol->size(); ++i){
          hs = &_hsCprCol->at(i);
          float  hs_t0 = hs->t0()._t0;
          if (std::abs(hs_t0 - ks_t0)<t0Toll){
            hasCprHelix = true;
            break;
          }
        }

        for (size_t i=0; i<_hsTprCol->size(); ++i){
          hs = &_hsTprCol->at(i);
          float  hs_t0 = hs->t0()._t0;
          if (std::abs(hs_t0 - ks_t0)<t0Toll){
            hasTprHelix = true;
            break;
          }
        }

        if ( hasCprHelix){
          _sumHist._hTrigInfo[20]->Fill(_nPOT);
        }
        if ( hasTprHelix ){
          _sumHist._hTrigInfo[21]->Fill(_nPOT);
        }

        for (size_t j=0; j<_effBits.size(); ++j){
          if (TrigResults->accept(_effBits[j])){
            _sumHist._hTrigInfo[17]->Fill(_nPOT);
            break;
          }
        }//end loop over the good effBits
      }//end loop over the tracks
    }

    if (Steps != NULL){
      const mu2e::StepPointMC* step(0);
      for ( size_t i=0; i<Steps->size(); ++i ){
        step = &Steps->at(i);
        int id = step->volumeId();
        if ( (id == 13) || (id == 14)){//Tracker front face
          if ( !(step->simParticle()->fromGenerator()) )                    continue;
          if (step->simParticle()->genParticle()->generatorId().id() != 43) continue;

          int   pdg = step->simParticle()->pdgId();
          float p   = step->momentum().mag();

          if ( (pdg == PDGCode::e_minus) && (p>=pCEleCuts[0]) && (p<=pCEleCuts[1])){
            _sumHist._hTrigInfo[18]->Fill(_nPOT);
            break;
          }else if ( (pdg == PDGCode::e_plus) && (p>=pCPosCuts[0]) && (p<=pCPosCuts[1])){
            _sumHist._hTrigInfo[18]->Fill(_nPOT);
            break;
          }
        }
      }//end loop over the virtual-det hits
    }
  }

  //--------------------------------------------------------------------------------
  void ReadTriggerInfo::analyze(const art::Event& event) {

    //get the number of POT
    _nPOT  = -1.;
    art::Handle<ProtonBunchIntensity> evtWeightH;
    event.getByLabel(_evtWeightTag, evtWeightH);
    if (evtWeightH.isValid()){
      _nPOT  = (double)evtWeightH->intensity();
    }
    _sumHist._hTrigInfo[19]->Fill(_nPOT);

    art::Handle<StepPointMCCollection> vdH;
    event.getByLabel(_vdTag, vdH);
    const mu2e::StepPointMCCollection*vdSteps(0);
    if (vdH.isValid()){
      vdSteps = vdH.product();
    }

    //get the TriggerResult
    std::ostringstream oss;
    oss << "TriggerResults::"<<_processName;
    art::InputTag const tag{oss.str()};
    auto const trigResultsH   = event.getValidHandle<art::TriggerResults>(tag);
    const art::TriggerResults*trigResults = trigResultsH.product();
    TriggerResultsNavigator   trigNavig(trigResults);

    //fill the histogram with the trigger bits
    //    for (unsigned i=0; i<trigResults->size(); ++i){
    for (unsigned int i=0; i< trigNavig.getTrigPaths().size(); ++i){
      //      if (trigResults->accept(i)){
      std::string path   = trigNavig.getTrigPathName(i);
      if(trigNavig.accepted(path)){
        for (size_t j=0; j<_trigPaths.size(); ++j){
          if (_trigPaths[j] == path){
            _sumHist._hTrigBits->Fill(j);
            break;
          }
        }
      }
    }

    for (unsigned int i=0; i< trigNavig.getTrigPaths().size(); ++i){
      std::string path   = trigNavig.getTrigPathName(i);
      size_t      pathID = trigNavig.findTrigPathID(path);
      if (trigNavig.accepted(path)) _sumHist._hTrigInfo[15]->Fill(pathID);
    }

    art::Handle<mu2e::HelixSeedCollection>  hsCprH;
    event.getByLabel(_hsCprTag, hsCprH);
    if (hsCprH.isValid()){
      _hsCprCol = hsCprH.product();
    }else {
      _hsCprCol = NULL;
    }

    art::Handle<mu2e::HelixSeedCollection>  hsTprH;
    event.getByLabel(_hsTprTag, hsTprH);
    if (hsTprH.isValid()){
      _hsTprCol = hsTprH.product();
    }else {
      _hsTprCol = NULL;
    }

    art::Handle<mu2e::KalSeedCollection>  ksH;
    const mu2e::KalSeedCollection*        ksCol(0);
    event.getByLabel(_ksTag, ksH);
    if (ksH.isValid()){
      ksCol = ksH.product();
    }
    art::Handle<mu2e::TrkQualCollection> trkQualH;
    const mu2e::TrkQualCollection*       trkQualCol(0);
    event.getByLabel(_trkQualTag, trkQualH);
    if (trkQualH.isValid()) {
      trkQualCol = trkQualH.product();
    }

    fillTrackEfficiencyHist(ksCol, trkQualCol, trigResults, vdSteps);

    //get the strawDigiMC truth if present
    art::Handle<mu2e::StrawDigiMCCollection> mcdH;
    event.getByLabel(_sdMCTag, mcdH);
    if (mcdH.isValid()) {
      _mcdigis = mcdH.product();
      _event   = &event;
    }else {
      _mcdigis = NULL;
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
      _chcol = NULL;
    }

    //get the CaloDigi Collection
    art::Handle<mu2e::CaloDigiCollection> cdH;
    event.getByLabel(_cdTag, cdH);
    const CaloDigiCollection* cdCol(0);
    if (cdH.isValid()) {
      cdCol = cdH.product();
    }

    //fill the general occupancy histogram
    fillOccupancyInfo   (_nTrackTrig+_nCaloTrig, sdCol, cdCol, _occupancyHist);

    std::vector<int>   trigFlagAll_index, trigFlag_index;

    art::Handle<TriggerInfo> hTrigInfoH;
    const mu2e::TriggerInfo* trigInfo(0);

    for (unsigned int i=0; i< _trigPaths.size(); ++i){
      string&path = _trigPaths.at(i);
      if (trigNavig.accepted(path)) {
        std::vector<std::string>      moduleNames = trigNavig.triggerModules(path);

        if(_diagLevel>0){
          printf("[ReadTriggerInfo::analyze] moduleNames size = %lu\n", moduleNames.size());
        }

        for (size_t j=0; j<moduleNames.size(); ++j){
          std::string  moduleLabel = moduleNames[j];
          if(_diagLevel>0){
            if (j==0) printf("[ReadTriggerInfo::analyze]      name      \n");
            printf("[ReadTriggerInfo::analyze] %10s\n", moduleLabel.c_str());
          }
          int          index_all(i);//0);
          int          index(i);//0);
          bool         passed(false);
          size_t       nTrigObj(0);
          //fill the Global Trigger bits info
          // findTrigIndex(_trigAll, moduleLabel, index_all);
          _trigAll[index_all].label  = moduleLabel;

          event.getByLabel(moduleLabel, hTrigInfoH);
          if (hTrigInfoH.isValid()){
            trigInfo = hTrigInfoH.product();
          }
          if ( moduleLabel.find(std::string("HSFilter")) != std::string::npos) {
            //            findTrigIndex(_trigHelix, moduleLabel, index);
            _trigHelix[index].label  = moduleLabel;
            _trigHelix[index].counts = _trigHelix[index].counts + 1;
            passed = true;
            nTrigObj=0;
            for (auto const hseed: trigInfo->helixes()){
              if(hseed){
                ++nTrigObj;
                fillHelixTrigInfo(index, hseed.get(), _helHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo(_nTrackTrig+index, sdCol, cdCol, _occupancyHist);
                }
              }
            }//end loop over the helix-collection
            _helHist._hHelInfo[i][120]->Fill(nTrigObj);

          }else if ( moduleLabel.find("KSFilter") != std::string::npos){
            //            findTrigIndex(_trigTrack, moduleLabel, index);
            _trigTrack[index].label  = moduleLabel;
            _trigTrack[index].counts = _trigTrack[index].counts + 1;
            passed = true;
            nTrigObj=0;
            for (auto const kseed: trigInfo->tracks()){
              if(kseed){
                ++nTrigObj;
                fillTrackTrigInfo(index, kseed.get(), _trkHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo(index, sdCol, cdCol, _occupancyHist);
                }
              }
            }//end loop over the kaseed-collection
            _trkHist._hTrkInfo[i][40]->Fill(nTrigObj);
            trigFlag_index.push_back(index_all);

          }else if ( moduleLabel.find("EventPrescale") != std::string::npos){
            //            findTrigIndex(_trigEvtPS, moduleLabel, index);
            _trigEvtPS[index].label  = moduleLabel;
            _trigEvtPS[index].counts = _trigEvtPS[index].counts + 1;
          }else if ( moduleLabel.find("caloCalibCosmic") != std::string::npos){
            //            findTrigIndex(_trigCaloCalib, moduleLabel, index);
            _trigCaloCalib[index].label  = moduleLabel;
            _trigCaloCalib[index].counts = _trigCaloCalib[index].counts + 1;
            passed = false;
            nTrigObj=0;
            for (auto const cluster : trigInfo->caloClusters()){
              if(cluster){
                ++nTrigObj;
                fillCaloCalibTrigInfo(index, cluster.get(), _caloCalibHist);
              }
            }//end loop over the cluster-collection
            trigFlag_index.push_back(index_all);
          }else if ( ( moduleLabel.find("caloMVANNCEFilter") != std::string::npos) || ( moduleLabel.find("caloPhotonFilter") != std::string::npos)){ //( (moduleLabel.find("caloMVACEFilter") != std::string::npos) || (moduleLabel.find("caloLHCEFilter") != std::string::npos) ){
            //            findTrigIndex(_trigCaloOnly, moduleLabel, index);
            _trigCaloOnly[index].label  = moduleLabel;
            _trigCaloOnly[index].counts = _trigCaloOnly[index].counts + 1;
            passed = true;
            nTrigObj=0;
            for (auto const clseed: trigInfo->caloClusters()){
              if(clseed){
                ++nTrigObj;
                fillCaloCalibTrigInfo(index, clseed.get(), _caloTSeedHist);
                if (passed) {
                  passed = false;
                  fillOccupancyInfo   (_nTrackTrig*2+index, sdCol, cdCol, _occupancyHist);
                }
              }
            }//end loop
            //_caloTSeedHist._hCaloOnlyInfo[i][20]->Fill(nTrigObj);
            trigFlag_index.push_back(index_all);
          }

          bool isCosmicHelix = (moduleLabel.find("HSFilter")!= std::string::npos) && (moduleLabel.find("CosmicHelix")!= std::string::npos);
          if ( (moduleLabel.find("caloPhotonFilter")!= std::string::npos) ||
               (moduleLabel.find("caloMVANNCEFilter")!= std::string::npos) ||
               (moduleLabel.find("TSFilter")       != std::string::npos) ||
               isCosmicHelix ){
            //            findTrigIndex(_trigFinal, moduleLabel, index);
            _trigFinal[index].label    = moduleLabel;
            _trigFinal[index].counts   = _trigFinal[index].counts + 1;
            _trigAll[index_all].counts = _trigAll[index_all].counts + 1;
            trigFlagAll_index.push_back(index_all);
          }
        }//end loop over the modules in a given trigger path
      }
    }

    //now fill the correlation matrix
    for (size_t i=0; i<trigFlagAll_index.size(); ++i){
      for (size_t j=0; j<trigFlagAll_index.size(); ++j){
        _sumHist._h2DTrigInfo[0]->Fill(trigFlagAll_index.at(i), trigFlagAll_index.at(j));
      }
    }

    if (trigFlagAll_index.size() == 1) _sumHist._hTrigInfo[10]->Fill(trigFlagAll_index.at(0));

  }

  void   ReadTriggerInfo::findTrigIndex(std::vector<trigInfo_> &Vec, std::string& ModuleLabel, int &Index){
    //reset the index value
    Index  = 0;
    size_t offset = 0;
    if ( ModuleLabel.find(std::string("tpr")) != std::string::npos) {
      offset = 10;
      Index  = 10;
    }
    for (size_t i=offset; i<Vec.size(); ++i){
      if ( (Vec[i].label == ModuleLabel) || (Vec[i].label == "") ) {
        Index = i;
        break;
      }// else if (Vec[i].label != ""){
      //         Index = i+1;
      // }
    }
  }

  void   ReadTriggerInfo::fillTrackTrigInfo(int TrkTrigIndex, const KalSeed*KSeed, trackInfoHist_   &Hist){
    GlobalConstantsHandle<ParticleDataList> pdt;
    //HelixTool helTool(KSeed->helix().get(), _tracker);

    int                nsh = (int)KSeed->hits().size();
    KalSegment const& fseg = KSeed->segments().front();

    double     ndof  = std::max(1.0,nsh - 5.0);
    double     p     = fseg.mom();
    double     chi2d = KSeed->chisquared()/ndof;
    double     pt    = p*std::cos(std::atan(fseg.helix().tanDip()));
    double     d0    = fseg.helix().d0();
    double     clE(-1.);
    if (KSeed->caloCluster()) clE = KSeed->caloCluster()->energyDep();
    //    double     nLoops    = helTool.nLoops();

    Hist._hTrkInfo[TrkTrigIndex][0]->Fill(p);
    Hist._hTrkInfo[TrkTrigIndex][1]->Fill(pt);
    Hist._hTrkInfo[TrkTrigIndex][2]->Fill(nsh);
    Hist._hTrkInfo[TrkTrigIndex][3]->Fill(d0);
    Hist._hTrkInfo[TrkTrigIndex][4]->Fill(chi2d);
    Hist._hTrkInfo[TrkTrigIndex][5]->Fill(clE);
    //    Hist._hTrkInfo[TrkTrigIndex][6]->Fill(nLoops);

    //add the MC info if available
    if (_mcdigis) {
      const mu2e::ComboHit*    hit(0), *hit_0(0);
      hit_0     = &_chcol->at(0);

      int                      loc(-1);
      std::vector<int>         hits_simp_id, hits_simp_index, hits_simp_z;

      for (int j=0; j<nsh; ++j){
        int  hitIndex  = int(KSeed->hits().at(j).index());
        hit            = &_chcol->at(hitIndex);
        loc            = hit - hit_0;
        const mu2e::StrawGasStep* step(0);
        if (loc > (int)_mcdigis->size()) {
          printf("[READTRIGGERINFO::fillTrackTrigInfo] loc = %d but MCDgis_size = %ld\n", loc, _mcdigis->size());
          continue;
        }
        const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(loc);
        if (sdmc->wireEndTime(mu2e::StrawEnd::cal) < sdmc->wireEndTime(mu2e::StrawEnd::hv)) {
          step = sdmc->strawGasStep(mu2e::StrawEnd::cal).get();
        }
        else {
          step = sdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
        }

        if (step) {
          art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();
          int sim_id        = simptr->id().asInt();

          hits_simp_id.push_back   (sim_id);
          hits_simp_index.push_back(loc);
          hits_simp_z.push_back(step->position().z());
        }
      }//end loop over the hits

      int     max(0), mostvalueindex(-1);//, mostvalue= hits_simp_id[0];
      float   dz_most(1e4);
      for (int k=0; k<(int)hits_simp_id.size(); ++k){
        int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
        if ( (co>0) &&  (co>max)) {
          float  dz      = std::fabs(hits_simp_z[k]);
          if (dz < dz_most){
            dz_most        = dz;
            max            = co;
            //            mostvalue      = hits_simp_id[k];
            mostvalueindex = hits_simp_index[k];
          }
        }
      }

      //finally, get the info of the first StrawDigi
      const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(mostvalueindex);
        art::Ptr<mu2e::SimParticle> const& simptr = sdmc->earlyStrawGasStep()->simParticle();
        int     pdg   = simptr->pdgId();
        art::Ptr<mu2e::SimParticle> mother = simptr;

        while(mother->hasParent())
          mother = mother->parent();
        //sim = mother.operator->();
        int      pdgM   = mother->pdgId();
        double   pXMC   = simptr->startMomentum().x();
        double   pYMC   = simptr->startMomentum().y();
        double   pZMC   = simptr->startMomentum().z();
        double   mass(-1.);//  = part->Mass();
        double   energy(-1.);// = sqrt(px*px+py*py+pz*pz+mass*mass);
        mass   = pdt->particle(pdg).mass();
        energy = sqrt(pXMC*pXMC+pYMC*pYMC+pZMC*pZMC+mass*mass);

        double   pTMC   = sqrt(pXMC*pXMC + pYMC*pYMC);
        double   pMC    = sqrt(pZMC*pZMC + pTMC*pTMC);

        CLHEP::Hep3Vector sp = simptr->startPosition();
        XYZVectorF origin;
        origin.SetX(sp.x()+3904);
        origin.SetY(sp.y());
        origin.SetZ(sp.z());
        double origin_r = sqrt(origin.x()*origin.x() + origin.y()*origin.y());
        double pz     = sqrt(p*p - pt*pt);

        //now fill the MC histograms
        Hist._hTrkInfo[TrkTrigIndex][10]->Fill(pMC);
        Hist._hTrkInfo[TrkTrigIndex][11]->Fill(pTMC);
        Hist._hTrkInfo[TrkTrigIndex][12]->Fill(pZMC);
        Hist._hTrkInfo[TrkTrigIndex][13]->Fill(p - pMC);
        Hist._hTrkInfo[TrkTrigIndex][14]->Fill(pt - pTMC);
        Hist._hTrkInfo[TrkTrigIndex][15]->Fill(pz - pZMC);
        Hist._hTrkInfo[TrkTrigIndex][16]->Fill(pdg);
        Hist._hTrkInfo[TrkTrigIndex][17]->Fill(origin.z());
        Hist._hTrkInfo[TrkTrigIndex][18]->Fill(origin_r);
        Hist._hTrkInfo[TrkTrigIndex][19]->Fill(pdgM);
        Hist._hTrkInfo[TrkTrigIndex][20]->Fill(energy);
    }
  }

  void     ReadTriggerInfo::fillHelixTrigInfoAdd     (int HelTrigIndex  , int MCMotherIndex, const HelixSeed* HSeed, helixInfoHist_         &Hist, MCInfo &TMPMCInfo){
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 0]->Fill(TMPMCInfo.pMC);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 1]->Fill(TMPMCInfo.p);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 2]->Fill(TMPMCInfo.d0);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 3]->Fill(TMPMCInfo.dpMC);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 4]->Fill(TMPMCInfo.dpTMC);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 5]->Fill(TMPMCInfo.dpZMC);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 6]->Fill(TMPMCInfo.pdg);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 7]->Fill(TMPMCInfo.origZ);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 8]->Fill(TMPMCInfo.origR);
    Hist._hHelInfo[HelTrigIndex][MCMotherIndex + 9]->Fill(TMPMCInfo.lambda);
  }


  void   ReadTriggerInfo::fillHelixTrigInfo(int HelTrigIndex, const HelixSeed*HSeed, helixInfoHist_  &Hist){
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

     //add the MC info if available
    if (_mcdigis) {
      //      const mu2e::ComboHit*    hit(0);
      std::vector<int>         hits_simp_id, hits_simp_index, hits_simp_z;
      float   minP(30.);
      for (int j=0; j<nch; ++j){
        std::vector<StrawDigiIndex> shids;
        HSeed->hits().fillStrawDigiIndices(j,shids);
        //        hit            = &HSeed->hits().at(j);

        for (size_t k=0; k<shids.size(); ++k) {
          const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(shids[k]);
          auto const& spmcp = sdmc->earlyStrawGasStep();
            art::Ptr<mu2e::SimParticle> const& simptr = spmcp->simParticle();
            int sim_id        = simptr->id().asInt();
            float   dz        = spmcp->position().z();// - trackerZ0;
            float   pMC       = std::sqrt(spmcp->momentum().mag2());
            if (pMC<minP)     continue;
            hits_simp_id.push_back   (sim_id);
            hits_simp_index.push_back(shids[k]);
            hits_simp_z.push_back(dz);
            break;
        }
      }//end loop over the hits

      if (hits_simp_id.size() == 0) return;

      int     max(0), mostvalueindex(-1);//, mostvalue= hits_simp_id[0];
      float   dz_most(1e4);
      for (int k=0; k<(int)hits_simp_id.size(); ++k){
        int co = (int)std::count(hits_simp_id.begin(), hits_simp_id.end(), hits_simp_id[k]);
        if ( (co>0) &&  (co>max)) {
          float  dz      = std::fabs(hits_simp_z[k]);
          if (dz < dz_most){
            dz_most        = dz;
            max            = co;
            //            mostvalue      = hits_simp_id[k];
            mostvalueindex = hits_simp_index[k];
          }
        }
      }

      //finally, get the info of the first StrawDigi
      const mu2e::StrawDigiMC* sdmc = &_mcdigis->at(mostvalueindex);
      auto const& spmcp = sdmc->earlyStrawGasStep();
        art::Ptr<mu2e::SimParticle> const& simptr = spmcp->simParticle();
        int     pdg   = simptr->pdgId();
        art::Ptr<mu2e::SimParticle> mother = simptr;

        while(mother->hasParent()) mother = mother->parent();
        int      pdgM   = mother->pdgId();
        double   pXMC   = spmcp->momentum().x();
        double   pYMC   = spmcp->momentum().y();
        double   pZMC   = spmcp->momentum().z();
        // double   mass(-1.);//  = part->Mass();
        // double   energy(-1.);// = sqrt(px*px+py*py+pz*pz+mass*mass);
        // mass   = pdt->particle(pdg).mass();
        // energy = sqrt(pXMC*pXMC+pYMC*pYMC+pZMC*pZMC+mass*mass);

        //need to check the mother of the particle
        //the possible cases we are interested are:
        // 1) IPA negative muon:
        // 2) atmospheric negative muon:
        // 3) atmospheric positive muon:
        // 4) photon
        // 5) proton
        // 6) neutron

        int   indexMother(-1);

        if (pdgM == PDGCode::mu_minus){ //negative muon
          XYZVectorF  mother_origin;
          mother_origin.SetX(mother->startPosition().x()+3904);
          mother_origin.SetY(mother->startPosition().y());
          mother_origin.SetZ(mother->startPosition().z());
          double mother_origin_r =  sqrt(mother_origin.x()*mother_origin.x() + mother_origin.y()*mother_origin.y());
          double IPA_z(7400.), IPA_z_tolerance(600.);
          // case 1: origin in the IPA and PDG-Id
          if ( (fabs(mother_origin_r) <= 300.5) &&
               (fabs(mother_origin_r) >= 299.5) &&
               (fabs(mother_origin.z() - IPA_z) < IPA_z_tolerance) &&
               (mother->startMomentum().px() < 1e-10) &&
               (mother->startMomentum().py() < 1e-10) &&
               (mother->startMomentum().pz() < 1e-10)){
            indexMother = 40;
          }else{//case 2: negative cosmic muon
            indexMother = 20;
          }
        }else if (pdgM == PDGCode::mu_plus){//case 3: positive muon
          indexMother = 30;
        }else if (pdgM == PDGCode::gamma){//case 4: photon
          indexMother = 50;
        }else if (pdgM == PDGCode::proton){//case 5: proton
          indexMother = 60;
        }else if (pdgM == PDGCode::n0){ //case 6: neutron
          indexMother = 70;
        }else if (pdgM == PDGCode::pi_minus){ //case 7: pi minus
          indexMother = 80;
        }else if (pdgM == PDGCode::pi_plus){ //case 8: pi plus
          indexMother = 90;
        }else if (pdgM == PDGCode::e_minus){ //case 9: electrons
          indexMother = 100;
        }else if (pdgM == PDGCode::e_plus){ //case 10: positrons
          indexMother = 110;
        }

        double   pTMC   = sqrt(pXMC*pXMC + pYMC*pYMC);
        double   pMC    = sqrt(pZMC*pZMC + pTMC*pTMC);

        CLHEP::Hep3Vector sp = simptr->startPosition();
        XYZVectorF origin;
        origin.SetX(sp.x()+3904);
        origin.SetY(sp.y());
        origin.SetZ(sp.z());
        double origin_r = sqrt(origin.x()*origin.x() + origin.y()*origin.y());
        // trackSeed->fOrigin1.SetXYZT(sp->x(),sp->y(),sp->z(),simptr->startGlobalTime());
        double pz     = sqrt(p*p - pt*pt);

        //now fill the MC histograms
        Hist._hHelInfo[HelTrigIndex][10]->Fill(pMC);
        Hist._hHelInfo[HelTrigIndex][11]->Fill(pTMC);
        Hist._hHelInfo[HelTrigIndex][12]->Fill(pZMC);
        Hist._hHelInfo[HelTrigIndex][13]->Fill(p - pMC);
        Hist._hHelInfo[HelTrigIndex][14]->Fill(pt - pTMC);
        Hist._hHelInfo[HelTrigIndex][15]->Fill(pz - pZMC);
        Hist._hHelInfo[HelTrigIndex][16]->Fill(pdg);
        Hist._hHelInfo[HelTrigIndex][17]->Fill(origin.z());
        Hist._hHelInfo[HelTrigIndex][18]->Fill(origin_r);
        Hist._hHelInfo[HelTrigIndex][19]->Fill(pdgM);

        //fill the "add" info
        if (indexMother>0){
          MCInfo tmpMCInfo;
          tmpMCInfo.pMC   = (pMC);
          tmpMCInfo.pTMC  = (pTMC);
          tmpMCInfo.pZMC  = (pZMC);
          tmpMCInfo.dpMC  = (p - pMC);
          tmpMCInfo.dpTMC = (pt - pTMC);
          tmpMCInfo.dpZMC = (pz - pZMC);
          tmpMCInfo.pdg   = (pdg);
          tmpMCInfo.origZ = (origin.z());
          tmpMCInfo.origR = (origin_r);
          tmpMCInfo.pdgM  = (pdgM);
          tmpMCInfo.lambda  = lambda;
          tmpMCInfo.d0    = d0;
          tmpMCInfo.p     = p;


          fillHelixTrigInfoAdd(HelTrigIndex, indexMother, HSeed, Hist, tmpMCInfo);
        }

        // Hist._hHelInfo[HelTrigIndex][20]->Fill(energy);
      }

  }
  //--------------------------------------------------------------------------------

  void   ReadTriggerInfo::fillCaloCalibTrigInfo(int ClCalibIndex, const CaloCluster*HCl, caloCalibrationHist_   &Hist){
    int        clsize    = HCl->size();
    double     energy    = HCl->energyDep();

    Hist._hCaloCalibInfo[ClCalibIndex][0]->Fill(energy);
    Hist._hCaloCalibInfo[ClCalibIndex][1]->Fill(clsize);
  }
  //--------------------------------------------------------------------------------

  // void   ReadTriggerInfo::fillCaloTrigSeedInfo(int Index, const CaloTrigSeed*HCl, caloCalibrationHist_     &Hist){
  //   int        clsize    = HCl->size();
  //   double     energy    = HCl->energyDep();

  //   Hist._hCaloCalibInfo[Index][0]->Fill(energy);
  //   Hist._hCaloCalibInfo[Index][1]->Fill(clsize);
  //   // Hist._hCaloOnlyInfo[Index][0]->Fill(HCl->epeak());
  //   // Hist._hCaloOnlyInfo[Index][1]->Fill(HCl->ring1max());
  //   // Hist._hCaloOnlyInfo[Index][2]->Fill(HCl->ring1max2());
  // }
  //--------------------------------------------------------------------------------

  void   ReadTriggerInfo::fillOccupancyInfo(int Index         , const StrawDigiCollection*SDCol, const CaloDigiCollection*CDCol, occupancyHist_   &Hist){
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
