
//
//
//
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
#include "art/Utilities/make_tool.h"
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

//Conditions
// #include "Offline/ConditionsService/inc/AcceleratorParams.hh"
// #include "Offline/ConditionsService/inc/ConditionsHandle.hh"

//Dataproducts
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh" //added
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh" //added
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"


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
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"


//ROOT
#include "TH1F.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include <TROOT.h>
#include "TLine.h"
#include "TEllipse.h"

#include <cmath>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <iostream>

namespace mu2e {

  class NPOTAnalysis : public art::EDAnalyzer {

  public:

    enum {
      kNEventHistsSets = 1,
      kNTimeClusterHistsSets = 1,
      kNHelixSeedHistsSets = 1
    };

    struct EventHists {
      TH1D* nTimeClusters;
      TH1D* nCaloHits; //initialize my new histograms that are filled per event here
      TH1D* caloEnergy;
      TH1D* nCaphriHits;
      TH1D* nProtonTCs;
      TH1D* nTrackerHits;

      TH1D* nPOT;
      TH1D* nProtonsDeltaFinder; //Number of POT from deltafinder module

      TH1D* recoNPOT_caloE; //1d RecoNPOT from RecoNPOTMaker_module.cc in CalPatRec

      TH1D* filteredRecoNPOT; //protons that passed the filter
      //2d hists! with nPOT
      TH2D* nTimeClusters2d;
      TH2D* nCaloHits2d;
      TH2D* caloEnergy2d;
      TH2D* nCaphriHits2d;
      TH2D* nProtonTCs2d;
      TH2D* nTrackerHits2d;

      //2d hists! wiht nProtonsDeltaFinder

      TH2D* nTimeClusters2df;
      TH2D* nCaloHits2df;
      TH2D* caloEnergy2df;
      TH2D* nCaphriHits2df;
      TH2D* nProtonTCs2df;
      TH2D* nTrackerHits2df;

      //2d hists! variable v variable
      TH2D* nPOTvnProtonTracks; //nPOT v nProtonsDeltaFinder
      TH2D* caloEH; //caloEnergy v nCaloHits
      TH2D* cTCs; //nTimeClusters v nProtonTCs (time clusters vrs from intensitycalo time clusters)
      TH2D* hitsTC; //nTrackerHits v nCaloHits
      TH2D* hitsCAPHRIvCal; //nCaphriHits v nCalHits
      TH2D* hitsCAPHRIvTrker; //nCaphriHits v nTrackerHits

      // Residual histograms:
      //TH1D* resid_nTimeClusters;
      // TH1D* resid_nCaloHits;
      TH1D* resid_caloE_reco; //residual histogram for nPOT predicted from CaloEnergy - nPOT from MCDataProducts divided by nPOT from MC
      //TH1D* extended_resid_caloE_reco;
      // TH1D* resid_nCaphriHits;
      // TH1D* resid_nProtonTCs;
      // TH1D* resid_nTrackerHits;

      //nPOT v residual histograms:
      // TH2D* resid_nTimeClusters2d;
      // TH2D* resid_nCaloHits2d;
      TH2D* resid_caloE_reco2d; //nPOT versus resid_caloEnergy normalized
      TH2D* resid_caloE_reco_v_caloE2d; // resid_caloEnergy versus caloE
      // TH2D* extended_resid_caloE_reco2d;
      // TH2D* resid_nCaphriHits2d;
      // TH2D* resid_nProtonTCs2d;
      // TH2D* resid_nTrackerHits2d;

    };

    struct TimeClusterHists {
      TH1D* dT;
      TH1D* t0;
    };

    struct HelixSeedHists{
      TH1D* radius;
    };

    struct Hists {
      EventHists* _EventHists[kNEventHistsSets]; //create an arry of  pointers of type EventHists called  _EventHists of size kNeventHistsSets (which here is 1)
      TimeClusterHists* _TimeClusterHists[kNTimeClusterHistsSets];
      HelixSeedHists* _HelixSeedHists[kNHelixSeedHistsSets];
    };

    struct tcData {
      double dT;
      double t0;
    };

    struct eventData {
      int nPOT_;
      int nProtonsDeltaFinder_; //number of protons from deltafinder module
      int nClusters;
      int nCaloHits_; //added bc given in IntensityInfoCalo.hh
      int caloEnergy_;
      int nCaphriHits_;
      int nProtonTCs_; //add in from intensityinfotimecluster.hh
      int nTrackerHits_; //from intensityinfotrackerhits.hh

      int recoNPOT_caloE_; //from RecoNPOTMaker_module.cc ProtonBunchIntensity.hh

      // Add residuals here:
      // double resid_nTimeClusters_;
      // double resid_nCaloHits_;
      double resid_caloE_reco_;
      // double resid_nCaphriHits_;
      // double resid_nProtonTCs_;
      // double resid_nTrackerHits_;

      bool passed; //for my filter
    };

    struct hsData{
      double radius;
    };

    explicit NPOTAnalysis(fhicl::ParameterSet const& pset);
    virtual ~NPOTAnalysis();

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    // This is called for each event.
    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir);
    void bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir);
    void bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir);

    void fillEventHistograms(EventHists* Hist);
    void fillTimeClusterHistograms(TimeClusterHists* Hist, size_t loopIndex);
    void fillHelixSeedHistograms(HelixSeedHists* Hist, size_t loopIndex);

    void bookHistograms(art::ServiceHandle<art::TFileService>& Tfs);
    void fillHistograms();

    void computeTimeClusterData();
    void computeEventData();
    void computeHelixSeedData();



  private:

    // parameter set stuff
    std::string                        _chCollTag;
    std::string                        _tcCollTag;
    std::string                        _hsCollTag;
    art::InputTag                      _evtWeightTag;


    Hists                              _hist;
    const mu2e::ComboHitCollection*    _chColl;
    const mu2e::TimeClusterCollection* _tcColl;
    const mu2e::ComboHit*              _ch;
    const art::Event*                  _event;
    const mu2e::HelixSeedCollection*   _hsColl; //this is initializing a pointer _hsColl of type constant HelixSeedCollection (from mu2e library)

    std::vector<tcData>                _tcData;
    std::vector<hsData>                _hsData; //makes a vector with each entry of type hsData and calls it _hsData

    //for filter
    std::string _processName;
    std::string _filterName;

    eventData                          _eventData;

  };



  NPOTAnalysis::NPOTAnalysis(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _chCollTag      (pset.get<string>("chCollTag", "TTmakePH")),
    _tcCollTag      (pset.get<string>("tcCollTag", "TTTZClusterFinder")),
    _hsCollTag      (pset.get<string>("hsCollTag", "TTAprHelixFinder")),
    _evtWeightTag   (pset.get<art::InputTag>("protonBunchIntensity" , "PBISim")),
    _processName    (pset.get<std::string>("processName", " S5Stn")), //idk if i did this correct at all, lumiStream is the name of the sequence my filter is in in CalPatTrkReco 's prolog.fcl
    _filterName     (pset.get<std::string>("filterName", "lumiStream"))  //changed RecoNPOTFilter to lumiStream
  {

  }
  //-----------------------------------------------------------------------------
  NPOTAnalysis::~NPOTAnalysis() {}

  //-----------------------------------------------------------------------------
  void NPOTAnalysis::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    char folder_name[20];
    TH1::AddDirectory(0);

    //-----------------------------------------------------------------------------
    // book event histograms
    //-----------------------------------------------------------------------------
    for (int i=0; i<kNEventHistsSets; i++) {
      sprintf(folder_name,"evt_%i",i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._EventHists[i] = new EventHists;
      bookEventHistograms(_hist._EventHists[i],&tfdir);
    }

    //-----------------------------------------------------------------------------
    // book time cluster histograms
    //-----------------------------------------------------------------------------
    for (int i=0; i<kNTimeClusterHistsSets; i++) {
      sprintf(folder_name,"tcl_%i",i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._TimeClusterHists[i] = new TimeClusterHists;
      bookTimeClusterHistograms(_hist._TimeClusterHists[i],&tfdir);
    }


    //-----------------------------------------------------------------------------
    // book helix seed histograms
    //-----------------------------------------------------------------------------

    for (int i=0; i<kNHelixSeedHistsSets; i++) {
      sprintf(folder_name, "hsc_%i",i);
      art::TFileDirectory tfdir = Tfs->mkdir(folder_name);
      _hist._HelixSeedHists[i] = new HelixSeedHists;
      bookHelixSeedHistograms(_hist._HelixSeedHists[i],&tfdir);
    }
  }
  //-----------------------------------------------------------------------------
  void NPOTAnalysis::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {

    //1D histograms!!!
    Hist->nTimeClusters = Dir->make<TH1D>("nTimeClusters" , "Number of Time Clusters; nTC ", 100, 0.0, 30.0 );
    Hist->nCaloHits     = Dir->make<TH1D>("nCaloHits", "Number of Calorimeter Hits; nCaloHits", 100, 0.0, 1300);
    Hist->caloEnergy    = Dir->make<TH1D>("caloEnergy", "Calorimeter Energy; E[MeV]", 100, 0.0, 7500);
    Hist->nCaphriHits   = Dir->make<TH1D>("nCAPHRIHits", "Number of CAPHRI Hits; nCAPHRIHits", 100, 0.0, 20);
    Hist->nProtonTCs    = Dir->make<TH1D>("nProtonTCs", "Number of Proton Time Clusters from Intensity Information; nProtonTC from Int", 100, 0.0, 100);
    Hist->nTrackerHits  = Dir->make<TH1D>("nTrackerHits", "Number of Tracker Hits from Intensity Information; nTrkHits from Int", 100, 0.0, 7500);

    Hist->nPOT          = Dir->make<TH1D>("nPOT", "Number of Protons on Target; nPOT", 100, 0.0, 1e8);
    Hist->nProtonsDeltaFinder        = Dir->make<TH1D>("nProtonsDeltaFinder", "Number of Proton tracks from DeltaFinder Module; nProtons ", 100, 0.0, 100);

    Hist->recoNPOT_caloE = Dir->make<TH1D>("recoNPOT_caloE", "Reconstructed nPOT from Calorimeter Energy; RecoNPOT from CaloE", 100, 0.0, 1e8);

    //filtered
    Hist->filteredRecoNPOT = Dir->make<TH1D>("filteredRecoNPOT", "Reconstructed nPOT that passed Filter; Filtered ReconPOT", 100, 0.0, 1e8);
    //2d histograms!!! nPOT
    Hist->nTimeClusters2d = Dir->make<TH2D>("nTimeClusters2d" , "Number of Protons on Target v Number of Time Clusters; nTC; nPOT ", 100, 0.0, 100.0, 100, 0.0, 1e8);
    Hist->nCaloHits2d     = Dir->make<TH2D>("nCaloHits2d", "Number of Protons on Target v Number of Calorimeter Hits; nCaloHits; nPOT", 100, 0.0, 1300, 100, 0.0, 1e8);
    Hist->caloEnergy2d    = Dir->make<TH2D>("caloEnergy2d", "Number of Protons on Target v Calorimeter Energy; E[MeV]; nPOT", 100, 0.0, 7500, 100, 0.0, 1e8);
    Hist->nCaphriHits2d   = Dir->make<TH2D>("nCAPHRIHits2d", "Number of Protons on Target v Number of CAPHRI Hits; nCAPHRIHits; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);
    Hist->nProtonTCs2d    = Dir->make<TH2D>("nProtonTCs2d", "Number of Protons on Target v Number of Proton Time Clusters from Intensity Information; nProtonTC from Int; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);
    Hist->nTrackerHits2d  = Dir->make<TH2D>("nTrackerHits2d", "Number of Protons on Target v Number of Tracker Hits from Intensity Information; nTrkHits from Int; nPOT", 100, 0.0, 7500, 100, 0.0, 1e8);

    //2d histograms!!!! nProtonsDeltaFinder
    Hist->nTimeClusters2df = Dir->make<TH2D>("nTimeClusters2df" , "Number of Protons from Delta Finder v Number of Time Clusters; nTC; nProtons from DeltaFinder ", 100, 0.0, 100.0, 100, 0.0, 100);
    Hist->nCaloHits2df     = Dir->make<TH2D>("nCaloHits2df", "Number of Protons from Delta Finder v Number of Calorimeter Hits; nCaloHits; nProtons from DeltaFinder", 100, 0.0, 1300, 100, 0.0, 100);
    Hist->caloEnergy2df    = Dir->make<TH2D>("caloEnergy2df", "Number of Protons from Delta Finder v Calorimeter Energy; E[MeV]; nProtons from DeltaFinder", 100, 0.0, 7500, 100, 0.0, 100);
    Hist->nCaphriHits2df   = Dir->make<TH2D>("nCAPHRIHits2df", "Number of Protons from Delta Finder v Number of CAPHRI Hits; nCAPHRIHits; nProtons from DeltaFinder", 100, 0.0, 100, 100, 0.0, 100);
    Hist->nProtonTCs2df    = Dir->make<TH2D>("nProtonTCs2df", "Number of Protons from Delta Finder v Number of Proton Time Clusters from Intensity Information; nProtonTC from Int; nProtons from DeltaFinder", 100, 0.0, 100, 100, 0.0, 100);
    Hist->nTrackerHits2df  = Dir->make<TH2D>("nTrackerHits2df", "Number of Protons from Delta Finder v Number of Tracker Hits from Intensity Information; nTrkHits from Int; nProtons from DeltaFinder", 100, 0.0, 7500, 100, 0.0, 100);

    //2d Hists!! variable v variable
    Hist->nPOTvnProtonTracks            = Dir->make<TH2D>("nPOTvnProtonTracks", "Number of Protons v Number of Protons from Delta Finder; nProtons from DF; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);
    Hist->caloEH           = Dir->make<TH2D>("caloEH", "Calorimeter Energy vs the number of Calorimeter Hits; nCaloHits; E[MeV]", 100, 0.0, 1300, 100, 0.0, 7500);
    Hist->cTCs             = Dir->make<TH2D>("cTCs", "Number of Time Clusters vs Number of Time Clusters from Intensity Information; nProtonTCs from IntensityInfo; nTCs", 100, 0.0, 100, 100, 0.0, 100);
    Hist->hitsTC           = Dir->make<TH2D>("hitsTC", "Number of Tracker Hits v number of Calorimeter Hits; nCaloHits; nTrackerHits", 100.0, 0.0, 1300, 20, 0.0, 7500);
    Hist->hitsCAPHRIvCal   = Dir->make<TH2D>("hitsCAPHRIvCal", "Number of CAPHRI Hits v Number of Calorimeter Hits; nCaloHits; nCAPHRIHits", 100, 0.0, 1300, 100, 0.0, 100);
    Hist->hitsCAPHRIvTrker = Dir->make<TH2D>("hitsCAPHRIvTrker", "Number of CAPHRI Hits v Number of Tracker Hits; nTrackerHits; nCAPHRIHits", 100, 0.0, 7500, 100, 0.0, 100);

    // Residual histograms:1D
    // Hist->resid_nTimeClusters = Dir->make<TH1D>("resid_nTimeClusters", "Residual between nPOT and predicted nPOT using the nTimeClusters; nPOT-PredictednPOT", 100, -1e6, 1e6);
    //Hist->resid_nCaloHits     = Dir->make<TH1D>("resid_nCaloHits", "Residual between nPOT and predicted nPOT using the nCaloHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);
    Hist->resid_caloE_reco    = Dir->make<TH1D>("resid_caloE_reco", "Normalized Residual:  nPOT - reconstucted  nPOT from Calorimeter Energy; ( nPOT-recoNPOT)/nPOT", 100, -0.600, 0.600); //eventually -.5,5 i hope
    //Hist->extended_resid_caloE_reco = Dir->make<TH1D>("extended_resid_caloE_reco", "Normalized Residual:  nPOT - reconstucted  nPOT from Calorimeter Energy; ( nPOT-recoNPOT)/nPOT", 100, -1.40, 0.40); //same but x-axis expanded to the left so no underflow
    // Hist->extended_resid_caloE_reco = Dir->make<TH1D>("extended_resid_caloE_reco", "Normalized Residual:  nPOT - reconstucted  nPOT from Calorimeter Energy; ( nPOT-recoNPOT)/nPOT", 100, -2.00, 2.00); //shouldnt need now that fixed!
    // Hist->resid_nCaphriHits   = Dir->make<TH1D>("resid_nCAPHRIHits", "Residual between nPOT and predicted nPOT using the nCAPHRIHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);
    // Hist->resid_nProtonTCs    = Dir->make<TH1D>("resid_nProtonTCs", "Residual between nPOT and predicted nPOT using the nProtonTCs;  nPOT-PredictednPOT", 100, -1e8, 1e8);
    // Hist->resid_nTrackerHits  = Dir->make<TH1D>("resid_nTrackerHits", "Residual between nPOT and predicted nPOT using the nTrackerHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);


    //nPOT v residual :2d ( nPOT v actualPOT-predictedPOT)
    // Hist->resid_nTimeClusters2d = Dir->make<TH2D>("resid_nTimeClusters2d", "nPOT v Residual of nPOT from nTimeClusters; nPOT-PredictednPOT; nPOT", 100, -1e6, 1e6, 100, 0.0, 1e8);
    // Hist->resid_nCaloHits2d     = Dir->make<TH2D>("resid_nCaloHits2d", "nPOT v Residual of nPOT from nCaloHits; nPOT-PredictednPOT; ", 100, -1e7, 1e7, 100, 0.0, 1e8);
    Hist->resid_caloE_reco2d    = Dir->make<TH2D>("resid_caloE_reco2d", "Normalized Residual of nPOT from Calorimeter Energy v nPOT;nPOT; (nPOT-recoNPOT)/nPOT", 100, 0.0, 1e8, 100, -0.600, 0.600);
    Hist->resid_caloE_reco_v_caloE2d = Dir->make<TH2D>("resid_caloE_reco_v_caloE2d", "Normalized Residual of nPOT from Calorimeter Energy v Calorimeter Energy; CaloE [MeV]; (nPOT-recoNPOT)/nPOT", 100, 0.0, 1e4, 100, -0.600, 0.600); //try 1e4 instead of 7500? so see where fit fails bc only should work for 7500 really idk just try
    // Hist->extended_resid_caloE_reco2d = Dir->make<TH2D>("extended_resid_caloE_reco2d", "nPOT v Normalized Residual of nPOT from Calorimeter Energy; (nPOT-recoNPOT)/nPOT; nPOT", 100, -1.40, 0.40, 100, 0.0, 1e8); //before giani recommendations
    // Hist->extended_resid_caloE_reco2d = Dir->make<TH2D>("extended_resid_caloE_reco2d", "Normalized Residual of nPOT from Calorimeter Energy vnPOT; nPOT; (nPOT-recoNPOT)/nPOT", 100, 0.0, 1e8, 100, -2.00, 2.00); //shouldnt need now that fixed
    //Hist->resid_nCaphriHits2d   = Dir->make<TH2D>("resid_nCAPHRIHits2d", "nPT v Residual of nPOT from nCAPHRIHits; nPOT-PredictednPOT; nPOT", 100, -1e7, 1e7, 100, 0.0, 1e8);
    // Hist->resid_nProtonTCs2d    = Dir->make<TH2D>("resid_nProtonTCs2d", "nPOT v Residual of nPOT from nProtonTCs; nPOT-PredictednPOT; nPOT", 100, -1e8, 1e8, 100, 0.0, 1e8);
    // Hist->resid_nTrackerHits2d  = Dir->make<TH2D>("resid_nTrackerHits2d", "nPOT v Residual of nPOT from nTrackerHits; nPOT-PredictednPOT; nPOT", 100, -1e7, 1e7, 100, 0.0, 1e8);
  }


  //-----------------------------------------------------------------------------
  void NPOTAnalysis::bookTimeClusterHistograms(TimeClusterHists* Hist, art::TFileDirectory* Dir) {

    Hist->dT = Dir->make<TH1D>("dT" , "width of time clusters", 40, 0.0, 500.0 );
    Hist->t0 = Dir->make<TH1D>("t0" , "time of cluster", 200, 0.0, 1800.0 );

  }

  //The one im adding
  //-----------------------------------------------------------------------------
  void NPOTAnalysis::bookHelixSeedHistograms(HelixSeedHists* Hist, art::TFileDirectory* Dir) {
    Hist->radius = Dir->make<TH1D>("radius", "radius of Helix", 50, 0.0, 500);
  }

  //--------------------------------------------------------------------------------//
  void NPOTAnalysis::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    bookHistograms(tfs);
  }

  //--------------------------------------------------------------------------------//
  void NPOTAnalysis::endJob(){}

  //--------------------------------------------------------------------------------//
  void NPOTAnalysis::beginRun(const art::Run & run){}

  //--------------------------------------------------------------------------------//
  void NPOTAnalysis::endSubRun(const art::SubRun& sr){}

  //--------------------------------------------------------------------------------//
  void NPOTAnalysis::fillHistograms() {

    //-----------------------------------------------------------------------------
    // fill event histograms
    //-----------------------------------------------------------------------------

    fillEventHistograms(_hist._EventHists[0]);



    //-----------------------------------------------------------------------------
    // fill time cluster histograms
    //-----------------------------------------------------------------------------
    for (size_t i=0; i<_tcData.size(); i++) {
      fillTimeClusterHistograms(_hist._TimeClusterHists[0], i);
    }

    //-----------------------------------------------------------------------------
    // fill helix seed histograms
    //-----------------------------------------------------------------------------

    for (size_t i=0; i<_hsData.size(); i++) {
      fillHelixSeedHistograms(_hist._HelixSeedHists[0], i);
    }
  }

  //-----------------------------------------------------------------------------
  void NPOTAnalysis::fillEventHistograms(EventHists* Hist) {

    //1D histograms!!
    Hist->nTimeClusters->Fill(_eventData.nClusters);
    Hist->nCaloHits->Fill(_eventData.nCaloHits_);
    Hist->caloEnergy->Fill(_eventData.caloEnergy_);
    Hist->nCaphriHits->Fill(_eventData.nCaphriHits_);
    Hist->nProtonTCs->Fill(_eventData.nProtonTCs_);
    Hist->nTrackerHits->Fill(_eventData.nTrackerHits_);

    Hist->nPOT->Fill(_eventData.nPOT_);
    Hist->nProtonsDeltaFinder->Fill(_eventData.nProtonsDeltaFinder_);

    //reconstructed npot
    Hist->recoNPOT_caloE->Fill(_eventData.recoNPOT_caloE_);

    //filtered
    if (_eventData.passed) {Hist->filteredRecoNPOT->Fill(_eventData.recoNPOT_caloE_); } //if the event passed the filter, then add it to the histogram
    //2d histograms!!! y is nPOT
    Hist->nTimeClusters2d->Fill(_eventData.nClusters, _eventData.nPOT_);
    Hist->nCaloHits2d->Fill(_eventData.nCaloHits_, _eventData.nPOT_);
    Hist->caloEnergy2d->Fill(_eventData.caloEnergy_, _eventData.nPOT_);
    Hist->nCaphriHits2d->Fill(_eventData.nCaphriHits_, _eventData.nPOT_);
    Hist->nProtonTCs2d->Fill(_eventData.nProtonTCs_, _eventData.nPOT_);
    Hist->nTrackerHits2d->Fill(_eventData.nTrackerHits_, _eventData.nPOT_);

    //2d histograms!!! y is nProtonsDeltaFinder
    Hist->nTimeClusters2df->Fill(_eventData.nClusters, _eventData.nProtonsDeltaFinder_);
    Hist->nCaloHits2df->Fill(_eventData.nCaloHits_, _eventData.nProtonsDeltaFinder_);
    Hist->caloEnergy2df->Fill(_eventData.caloEnergy_, _eventData.nProtonsDeltaFinder_);
    Hist->nCaphriHits2df->Fill(_eventData.nCaphriHits_, _eventData.nProtonsDeltaFinder_);
    Hist->nProtonTCs2df->Fill(_eventData.nProtonTCs_, _eventData.nProtonsDeltaFinder_);
    Hist->nTrackerHits2df->Fill(_eventData.nTrackerHits_, _eventData.nProtonsDeltaFinder_);

    //2d Hists variable v variable
    Hist->nPOTvnProtonTracks->Fill(_eventData.nProtonsDeltaFinder_, _eventData.nPOT_);
    Hist->caloEH->Fill(_eventData.nCaloHits_, _eventData.caloEnergy_);
    Hist->cTCs->Fill(_eventData.nProtonTCs_, _eventData.nClusters);
    Hist->hitsTC->Fill(_eventData.nCaloHits_, _eventData.nTrackerHits_);
    Hist->hitsCAPHRIvCal->Fill(_eventData.nCaloHits_, _eventData.nCaphriHits_);
    Hist->hitsCAPHRIvTrker->Fill(_eventData.nTrackerHits_, _eventData.nCaphriHits_);


    // fill residual histograms
    // Hist->resid_nTimeClusters->Fill(_eventData.resid_nTimeClusters_);
    //Hist->resid_nCaloHits->Fill(_eventData.resid_nCaloHits_);
    Hist->resid_caloE_reco->Fill(_eventData.resid_caloE_reco_);
    // Hist->extended_resid_caloE_reco->Fill(_eventData.resid_caloE_reco_);
    //Hist->resid_nCaphriHits->Fill(_eventData.resid_nCaphriHits_);
    // Hist->resid_nProtonTCs->Fill(_eventData.resid_nProtonTCs_);
    // Hist->resid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_);

    //fill the 2d nPOT v residual histograms
    // Hist->resid_nTimeClusters->Fill(_eventData.resid_nTimeClusters_, _eventData.nPOT_);
    //Hist->resid_nCaloHits->Fill(_eventData.resid_nCaloHits_, _eventData.nPOT_);
    Hist->resid_caloE_reco2d->Fill( _eventData.nPOT_, _eventData.resid_caloE_reco_); //residuals on y-axis
    Hist->resid_caloE_reco_v_caloE2d->Fill(_eventData.caloEnergy_, _eventData.resid_caloE_reco_);
    //Hist->extended_resid_caloE_reco2d->Fill(_eventData.resid_caloE_reco_, _eventData.nPOT_);
    // Hist->extended_resid_caloE_reco2d->Fill(_eventData.nPOT_, _eventData.resid_caloE_reco_); //residual v nPOT ask giani advised.

    //Hist->resid_nCaphriHits->Fill(_eventData.resid_nCaphriHits_, _eventData.nPOT_);
    //Hist->resid_nProtonTCs->Fill(_eventData.resid_nProtonTCs_, _eventData.nPOT_);
    //Hist->resid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_, _eventData.nPOT_);


  }

  //-----------------------------------------------------------------------------
  void NPOTAnalysis::fillTimeClusterHistograms(TimeClusterHists* Hist, size_t loopIndex) {

    // fill dT plot
    double deltaTime = _tcData.at(loopIndex).dT;
    Hist->dT->Fill(deltaTime);
    // fill t0 plot
    double timeZero = _tcData.at(loopIndex).t0;
    Hist->t0->Fill(timeZero);

  }


  //-----------------------------------------------------------------------------
  void NPOTAnalysis::fillHelixSeedHistograms(HelixSeedHists* Hist, size_t loopIndex) {
    //fill radius plot
    double radiusdata = _hsData.at(loopIndex).radius;
    Hist->radius->Fill(radiusdata);
  }
  //--------------------------------------------------------------------------------
  void NPOTAnalysis::computeEventData() {

    // for (size_t i=0; i<_chCollph->size(); i++) {



    // compute event level stuff
    _eventData.nClusters = (int)_tcColl->size();

    //find the residual btw the predicted nPOT( from observed _eventData.observable plugged into  fit function defined) and the actual nPOT. residual = actual-pred
    /* if (_eventData.nClusters > 0) {
       double pred = nTimeClusters2d(static_cast<double>(_eventData.nClusters));
       _eventData.resid_nTimeClusters_ = static_cast<double>(_eventData.nPOT_) - pred;
       } else {
       _eventData.resid_nTimeClusters_ = 0.0;
       }

       if (_eventData.nCaloHits_ > 0) {
       double pred = nCaloHits2d(static_cast<double>(_eventData.nCaloHits_));
       _eventData.resid_nCaloHits_ = static_cast<double>(_eventData.nPOT_) - pred;
       } else {
       _eventData.resid_nCaloHits_ = 0.0;
       }
    */



    if (_eventData.nPOT_ > 0 && _eventData.recoNPOT_caloE_ > 0) {
      double residual = _eventData.nPOT_ - _eventData.recoNPOT_caloE_;
      _eventData.resid_caloE_reco_ = residual / _eventData.nPOT_;
    } else {
      _eventData.resid_caloE_reco_ = 0.0;
      std::cout << "Failed to calculate normalized residual: nPOT = "
                << _eventData.nPOT_ << ", recoNPOT = " << _eventData.recoNPOT_caloE_ << std::endl;
    }

  }
  /*
    if (_eventData.nCaphriHits_ > 0) {
    double pred = nCaphriHits2d(static_cast<double>(_eventData.nCaphriHits_));
    _eventData.resid_nCaphriHits_ = static_cast<double>(_eventData.nPOT_) - pred;
    } else {
    _eventData.resid_nCaphriHits_ = 0.0;
    }

    if (_eventData.nProtonTCs_ > 0) {
    double pred = nProtonTCs2d(static_cast<double>(_eventData.nProtonTCs_));
    _eventData.resid_nProtonTCs_ = static_cast<double>(_eventData.nPOT_) - pred;
    } else {
    _eventData.resid_nProtonTCs_ = 0.0;
    }

    if (_eventData.nTrackerHits_ > 0) {
    double pred = nTrackerHits2d(static_cast<double>(_eventData.nTrackerHits_));
    _eventData.resid_nTrackerHits_ = static_cast<double>(_eventData.nPOT_) - pred;
    } else {
    _eventData.resid_nTrackerHits_ = 0.0;

    }
    }

  */
  //--------------------------------------------------------------------------------
  void NPOTAnalysis::computeTimeClusterData() {




    // Get the TimeClusterCollection
    //auto tcCollHandle = _event->getValidHandle<TimeClusterCollection>(_chCollTag);
    //_tcColl = const_cast<TimeClusterCollection*>(tcCollHandle.product());


    // compute dT for each cluster
    for (size_t i=0; i<_tcColl->size(); i++) {
      tcData clusterData;
      clusterData.t0 = _tcColl->at(i)._t0._t0;
      double minT = 0.0;
      double maxT = 0.0;
      for (size_t j=0; j<_tcColl->at(i)._strawHitIdxs.size(); j++) {
        int chIndex = _tcColl->at(i)._strawHitIdxs[j];
        double chTime = _chColl->at(chIndex).correctedTime();
        if (j==0) {
          minT = chTime;
          maxT = chTime;
        }
        if (chTime<minT) {minT = chTime;}
        if (chTime>maxT) {maxT = chTime;}
      }
      clusterData.dT = maxT - minT;
      _tcData.push_back(clusterData);
    }

  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::computeHelixSeedData() {


    //loop through all of the helices and grab the data members of interest
    for (size_t i=0; i<_hsColl->size(); i++) {
      hsData radiusData;
      radiusData.radius = _hsColl->at(i).helix().radius();
      _hsData.push_back(radiusData);
    }

  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::analyze(const art::Event& event) {

    // get event
    _event = &event;


    //LOAD IN WHETHER PASSED FILTER OR NOT from RecoNPOTFilter_module.cc
    std::ostringstream oss;
    oss << "TriggerResults::" << _processName; //fcl parameter
    art::InputTag const tag{oss.str()};

    /*DEBUGGING TO MAKE SURE PROCESS AND FILTER NAME ARE VALID
    std::cout << "[DEBUG] Looking for TriggerResults with tag: '" << tag.encode() << "'" << std::endl;
    std::cout << "[DEBUG] Process name: '" << _processName << "'" << std::endl;
    std::cout << "[DEBUG] Filter name: '" << _filterName << "'" << std::endl;
    */
    try {
      auto const trigResultsH = event.getValidHandle<art::TriggerResults>(tag);
      if (trigResultsH.isValid()) {
        // std::cout << "[DEBUG] TriggerResults handle IS VALID" << std::endl;
        const art::TriggerResults* trigResults = trigResultsH.product();
        TriggerResultsNavigator trigNavig(trigResults);
        /*DEBUGGING TO SEE AVALIBLE PATH NAMES-----------------

        std::cout << "[DEBUG] Paths and acceptance from TriggerResultsNavigator:" << std::endl;
        for (size_t i = 0; i < trigResults->size(); ++i) {
          std::string const& pathName = trigNavig.getTrigPathName(i);
          bool accepted = trigNavig.accepted(pathName);
          std::cout << "  [" << i << "] Path: " << pathName << std::endl;
        }
        //---------------------------------------------------*/
        _eventData.passed = trigNavig.accepted(_filterName);

        /*DEBUGGING TO SEE IF PASSED FILTER
        std::cout << "[NPOTAnalysis] Filter check successful. Event " << event.id()
                  << " passed filter '" << _filterName << "': "
                  << (_eventData.passed ? "YES" : "NO") << std::endl;
        */
      } else {
        std::cout << "[NPOTAnalysis] TriggerResults handle is not valid!" << std::endl;
        _eventData.passed = false;
      }
    } catch (const std::exception& e) {
      std::cout << "[NPOTAnalysis] Could not retrieve TriggerResults with tag '"
                << tag.encode() << "'. Error: " << e.what() << std::endl;
      _eventData.passed = false;
    }




    //RECONSTRUCTED NPOT FOR ANALYSIS ON FIT AND FOR RecoNPOTMaker_module.cc
    art::Handle<ProtonBunchIntensity> recoNPOTH;
    _event->getByLabel("RecoNPOTMaker",recoNPOTH);
    if(recoNPOTH.isValid()) {
      _eventData.recoNPOT_caloE_ = recoNPOTH->intensity();
      // std::cout << "Event ID: " << event.id()
      //<< "  recoNPOT_caloE_ = " << _eventData.recoNPOT_caloE_ << std::endl;

    }
    else { std::cout << "[NPOTAnalysis] Could not retrieve RecoNPOT from RecoNPOTMaker!" << std::endl; }



    //get the EventData from ProtonBunchIntensity
    art::Handle<ProtonBunchIntensity> evtWeightH;
    _event->getByLabel(_evtWeightTag, evtWeightH);
    if (evtWeightH.isValid()) {_eventData.nPOT_  = evtWeightH->intensity();}

    //get the EventData from IntensityInfoCalo
    art::Handle<IntensityInfoCalo> caloInt;
    _event->getByLabel("CaloHitMakerFast", caloInt);
    if (caloInt.isValid()){
      _eventData.nCaloHits_   = caloInt->nCaloHits();
      _eventData.caloEnergy_  = caloInt->caloEnergy();
      _eventData.nCaphriHits_ = caloInt->nCaphriHits();
    }

    art::Handle<IntensityInfoTimeCluster> tcInt;
    _event->getByLabel("TTTZClusterFinder", tcInt);
    if (tcInt.isValid()) { _eventData.nProtonTCs_  = tcInt->nProtonTCs();}

    //Get the nProtonsDeltaFinder from IntensityInfoTimeCluster
    art::Handle<IntensityInfoTimeCluster> tcIntph;
    _event->getByLabel("TTflagPH", tcIntph);
    if (tcIntph.isValid()) { _eventData.nProtonsDeltaFinder_ = tcIntph->nProtonTCs();}
    else { std::cout << "tcIntph is not valid, not filling nProtonsDeltaFinder_" << std::endl;}

    art::Handle<IntensityInfoTrackerHits> thInt;
    _event->getByLabel("TTmakeSH", thInt);
    if (thInt.isValid()) { _eventData.nTrackerHits_ = thInt->nTrackerHits();}





    // get the ComboHitCollection
    art::Handle<mu2e::ComboHitCollection> chCollHandle;
    _event->getByLabel(_chCollTag, chCollHandle);
    if (chCollHandle.isValid())
      {
        _chColl = chCollHandle.product();
      }
    else {_chColl = NULL;}


    //get the TimeClusterCollection
    art::Handle<mu2e::TimeClusterCollection> tcCollHandle;
    _event->getByLabel(_tcCollTag, tcCollHandle);
    if (tcCollHandle.isValid()) {_tcColl = tcCollHandle.product();}
    else {_tcColl = NULL;}

    //get the HelixSeedCollection
    art::Handle<mu2e::HelixSeedCollection> hsCollHandle;
    _event->getByLabel(_hsCollTag, hsCollHandle);
    if (hsCollHandle.isValid()) {_hsColl = hsCollHandle.product();}
    else {_hsColl = NULL;}


    computeEventData();
    // compute data that will be plotted
    if(_tcColl != NULL){ computeTimeClusterData();}
    //radius idk man
    if(_hsColl != NULL){ computeHelixSeedData();}



    // fill histograms
    fillHistograms();

    _tcData.clear();
    _hsData.clear();

  }


}
DEFINE_ART_MODULE(mu2e::NPOTAnalysis)
