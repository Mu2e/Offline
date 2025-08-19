//
// Analyze MC and reco N(POT) info
// Original author: Hope Applegate, 2025
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Data products
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/RecoProtonBunchIntensity.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

//MC data products
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"

//Utilities
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"


//ROOT
#include "TH1.h"
#include "TH2.h"

//c++
#include <string>
#include <vector>
#include <iostream>

namespace mu2e {

  class NPOTAnalysis : public art::EDAnalyzer {

  public:

    enum {
      kNEventHistsSets = 1
    };

    struct EventHists {
      //-------------1D Histograms of all the Observables-------------
      TH1D* nTimeClusters;
      TH1D* nCaloHits;
      TH1D* caloEnergy;
      TH1D* nCaphriHits;
      TH1D* nProtonTCs;
      TH1D* nTrackerHits;
      TH1D* nProtonsDeltaFinder; //Number of Protons from DeltaFinder module

      TH1D* nPOT; //MC Truth  nPOT

      TH1D* recoNPOT_caloE; //Produced in RecoNPOTMaker_module.cc in CalPatRec

      TH1D* filteredRecoNPOT; //Passed RecoNPOTFilter_module.cc

      //-------------Observable v. NPOT-------------
      TH2D* nTimeClusters2d;
      TH2D* nCaloHits2d;
      TH2D* caloEnergy2d;
      TH2D* nCaphriHits2d;
      TH2D* nProtonTCs2d;
      TH2D* nTrackerHits2d;
      TH2D* nPOTvnProtonTracks; //nPOT v nProtonsDeltaFinder

      //-------------Observable v. nProtons from Delta Finder Module-------------
      TH2D* nTimeClusters2df;
      TH2D* nCaloHits2df;
      TH2D* caloEnergy2df;
      TH2D* nCaphriHits2df;
      TH2D* nProtonTCs2df;
      TH2D* nTrackerHits2df;

      //-------------Observable v Observable-------------
      TH2D* caloEH; //caloEnergy v nCaloHits
      TH2D* cTCs; //nTimeClusters v nProtonTCs (from IntensityCalo Time Clusters)
      TH2D* hitsTC; //nTrackerHits v nCaloHits
      TH2D* hitsCAPHRIvCal; //nCaphriHits v nCalHits
      TH2D* hitsCAPHRIvTrker; //nCaphriHits v nTrackerHits

      //--------------Residual histograms: (nPOT predicted from Observable-nPOT)/(nPOT)--------------
      //naming scheme: resid_Observable
      TH1D* resid_nTimeClusters;
      TH1D* resid_nCaloHits;
      TH1D* resid_caloE_reco;
      TH1D* resid_nCaphriHits;
      TH1D* resid_nProtonTCs;
      TH1D* resid_nTrackerHits;

      //--------------nPOT v residual histograms:--------------
      TH2D* resid_nTimeClusters2d;
      TH2D* resid_nCaloHits2d;
      TH2D* resid_caloE_reco2d; //nPOT versus resid_caloEnergy normalized
      TH2D* resid_nCaphriHits2d;
      TH2D* resid_nProtonTCs2d;
      TH2D* resid_nTrackerHits2d;

      //-------------Observable v Residual-------------
      TH2D* resid_caloE_reco_v_caloE2d; // resid_caloEnergy versus caloE

    };

    struct Hists {
      EventHists* _EventHists[kNEventHistsSets];
    };

    struct eventData {
      long long nPOT_;
      //-------------Observables-------------
      int nProtonsDeltaFinder_;
      int nClusters_;
      int nCaloHits_;
      int caloEnergy_;
      int nCaphriHits_;
      int nProtonTCs_;
      int nTrackerHits_;
      //-------------Made in RecoNPOTMaker_module.cc as <RecoProtonBunchIntensity>
      long long recoNPOT_caloE_;

      //-------------Residuals-------------
      double resid_nTimeClusters_;
      double resid_nCaloHits_;
      double resid_caloE_reco_;
      double resid_nCaphriHits_;
      double resid_nProtonTCs_;
      double resid_nTrackerHits_;

      //-------------Passed RecoNPOTFilter_module.cc-------------
      bool passed_;

      eventData() { reset(); }

      void reset() {
        nPOT_ = 0;
        nProtonsDeltaFinder_ = 0;
        nClusters_ = 0;
        nCaloHits_ = 0;
        caloEnergy_ = 0;
        nCaphriHits_ = 0;
        nProtonTCs_ = 0;
        nTrackerHits_ = 0;
        recoNPOT_caloE_ = 0;
        resid_nTimeClusters_ = 0;
        resid_nCaloHits_ = 0;
        resid_caloE_reco_ = 0;
        resid_nCaphriHits_ = 0;
        resid_nProtonTCs_ = 0;
        resid_nTrackerHits_ = 0;
        passed_ = 0;
      }
    };

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> caloTag   {Name("caloTag")    , Comment("Calo cluster tag")};
      fhicl::Atom<art::InputTag> thCollTag {Name("thCollTag")  , Comment("Straw hit collection tag")};
      fhicl::Atom<art::InputTag> dfCollTag {Name("dfCollTag")  , Comment("Delta finder tag")};
      fhicl::Atom<art::InputTag> tcCollTag {Name("tcCollTag")  , Comment("Time cluster collection")};
      fhicl::Atom<art::InputTag> recoPOTTag{Name("recoPOTTag") , Comment("Reco POT tag")};
      fhicl::Atom<art::InputTag> MCPOTTag  {Name("MCPOTTag")   , Comment("MC POT tag")};
      fhicl::Atom<std::string>   process   {Name("processName"), Comment("Name of this process")};
      fhicl::Atom<std::string>   filter    {Name("filterName") , Comment("Name of the N(POT) filter")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit NPOTAnalysis(const Parameters& conf);
    virtual ~NPOTAnalysis();

    virtual void beginJob();
    virtual void endJob();
    virtual void endSubRun(const art::SubRun& sr);

    virtual void analyze(const art::Event& e);
    virtual void beginRun(const art::Run & run);

    void bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir);
    void fillEventHistograms(EventHists* Hist);

    void bookHistograms(art::ServiceHandle<art::TFileService>& Tfs);
    void fillHistograms();

    void retrieveData();
    void computeEventData();



  private:

    art::InputTag                      _caloTag;
    art::InputTag                      _thCollTag;
    art::InputTag                      _dfCollTag;
    art::InputTag                      _tcCollTag;
    art::InputTag                      _recoPOTTag;
    art::InputTag                      _MCPOTTag;
    std::string                        _processName;
    std::string                        _filterName;

    const art::Event*                  _event;
    const mu2e::TimeClusterCollection* _tcColl;

    eventData                          _eventData;
    Hists                              _hist;

  };

  //-----------------------------------------------------------------------------
  NPOTAnalysis::NPOTAnalysis(const Parameters& conf) : art::EDAnalyzer{conf}
    , _caloTag       (conf().caloTag())
    , _thCollTag     (conf().thCollTag())
    , _dfCollTag     (conf().dfCollTag())
    , _tcCollTag     (conf().tcCollTag())
    , _recoPOTTag    (conf().recoPOTTag())
    , _MCPOTTag      (conf().MCPOTTag())
    , _processName   (conf().process())
    , _filterName    (conf().filter())
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
  }

  //-----------------------------------------------------------------------------
  void NPOTAnalysis::bookEventHistograms(EventHists* Hist, art::TFileDirectory* Dir) {

    //-------------1D Observables-------------
    Hist->nTimeClusters          = Dir->make<TH1D>("nTimeClusters" , "Number of Time Clusters; nTC ", 30, 0.0, 30.0 );
    Hist->nCaloHits              = Dir->make<TH1D>("nCaloHits", "Number of Calorimeter Hits; nCaloHits", 100, 0.0, 1300);
    Hist->caloEnergy             = Dir->make<TH1D>("caloEnergy", "Calorimeter Energy; E[MeV]", 100, 0.0, 7500);
    Hist->nCaphriHits            = Dir->make<TH1D>("nCAPHRIHits", "Number of CAPHRI Hits; nCAPHRIHits", 20, 0.0, 20);
    Hist->nProtonTCs             = Dir->make<TH1D>("nProtonTCs", "Number of Proton Time Clusters from Intensity Information; nProtonTC from Int", 100, 0.0, 100);
    Hist->nTrackerHits           = Dir->make<TH1D>("nTrackerHits", "Number of Tracker Hits from Intensity Information; nTrkHits from Int", 100, 0.0, 7500);
    Hist->nProtonsDeltaFinder    = Dir->make<TH1D>("nProtonsDeltaFinder", "Number of Proton tracks from DeltaFinder Module; nProtons ", 100, 0.0, 100);

    //-------------MC Truth nPOT-------------
    Hist->nPOT                   = Dir->make<TH1D>("nPOT", "Number of Protons on Target; nPOT", 100, 0.0, 1e8);

    //-------------Reconstructed nPOT from Observable-------------
    Hist->recoNPOT_caloE         = Dir->make<TH1D>("recoNPOT_caloE", "Reconstructed nPOT from Calorimeter Energy; RecoNPOT from CaloE", 100, 0.0, 1e8);

    //-------------Filtered-------------
    Hist->filteredRecoNPOT       = Dir->make<TH1D>("filteredRecoNPOT", "Reconstructed nPOT that passed Filter; Filtered ReconPOT", 100, 0.0, 1e8);

    //-------------Observable v nPOT-------------
    Hist->nTimeClusters2d        = Dir->make<TH2D>("nTimeClusters2d" , "Number of Protons on Target v Number of Time Clusters; nTC; nPOT ", 100, 0.0, 100.0, 100, 0.0, 1e8);
    Hist->nCaloHits2d            = Dir->make<TH2D>("nCaloHits2d", "Number of Protons on Target v Number of Calorimeter Hits; nCaloHits; nPOT", 100, 0.0, 1300, 100, 0.0, 1e8);
    Hist->caloEnergy2d           = Dir->make<TH2D>("caloEnergy2d", "Number of Protons on Target v Calorimeter Energy; E[MeV]; nPOT", 100, 0.0, 7500, 100, 0.0, 1e8);
    Hist->nCaphriHits2d          = Dir->make<TH2D>("nCAPHRIHits2d", "Number of Protons on Target v Number of CAPHRI Hits; nCAPHRIHits; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);
    Hist->nProtonTCs2d           = Dir->make<TH2D>("nProtonTCs2d", "Number of Protons on Target v Number of Proton Time Clusters from Intensity Information; nProtonTC from Int; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);
    Hist->nTrackerHits2d         = Dir->make<TH2D>("nTrackerHits2d", "Number of Protons on Target v Number of Tracker Hits from Intensity Information; nTrkHits from Int; nPOT", 100, 0.0, 7500, 100, 0.0, 1e8);
    Hist->nPOTvnProtonTracks     = Dir->make<TH2D>("nPOTvnProtonTracks", "Number of Protons v Number of Protons from Delta Finder; nProtons from DF; nPOT", 100, 0.0, 100, 100, 0.0, 1e8);

    //-------------Observable v nProtons from DeltaFinder-------------
    Hist->nTimeClusters2df       = Dir->make<TH2D>("nTimeClusters2df" , "Number of Protons from Delta Finder v Number of Time Clusters; nTC; nProtons from DeltaFinder ", 100, 0.0, 100.0, 100, 0.0, 100);
    Hist->nCaloHits2df           = Dir->make<TH2D>("nCaloHits2df", "Number of Protons from Delta Finder v Number of Calorimeter Hits; nCaloHits; nProtons from DeltaFinder", 100, 0.0, 1300, 100, 0.0, 100);
    Hist->caloEnergy2df          = Dir->make<TH2D>("caloEnergy2df", "Number of Protons from Delta Finder v Calorimeter Energy; E[MeV]; nProtons from DeltaFinder", 100, 0.0, 7500, 100, 0.0, 100);
    Hist->nCaphriHits2df         = Dir->make<TH2D>("nCAPHRIHits2df", "Number of Protons from Delta Finder v Number of CAPHRI Hits; nCAPHRIHits; nProtons from DeltaFinder", 100, 0.0, 100, 100, 0.0, 100);
    Hist->nProtonTCs2df          = Dir->make<TH2D>("nProtonTCs2df", "Number of Protons from Delta Finder v Number of Proton Time Clusters from Intensity Information; nProtonTC from Int; nProtons from DeltaFinder", 100, 0.0, 100, 100, 0.0, 100);
    Hist->nTrackerHits2df        = Dir->make<TH2D>("nTrackerHits2df", "Number of Protons from Delta Finder v Number of Tracker Hits from Intensity Information; nTrkHits from Int; nProtons from DeltaFinder", 100, 0.0, 7500, 100, 0.0, 100);

    //-------------Observable v Observable-------------

    Hist->caloEH                 = Dir->make<TH2D>("caloEH", "Calorimeter Energy vs the number of Calorimeter Hits; nCaloHits; E[MeV]", 100, 0.0, 1300, 100, 0.0, 7500);
    Hist->cTCs                   = Dir->make<TH2D>("cTCs", "Number of Time Clusters vs Number of Time Clusters from Intensity Information; nProtonTCs from IntensityInfo; nTCs", 100, 0.0, 100, 100, 0.0, 100);
    Hist->hitsTC                 = Dir->make<TH2D>("hitsTC", "Number of Tracker Hits v number of Calorimeter Hits; nCaloHits; nTrackerHits", 100.0, 0.0, 1300, 20, 0.0, 7500);
    Hist->hitsCAPHRIvCal         = Dir->make<TH2D>("hitsCAPHRIvCal", "Number of CAPHRI Hits v Number of Calorimeter Hits; nCaloHits; nCAPHRIHits", 100, 0.0, 1300, 100, 0.0, 100);
    Hist->hitsCAPHRIvTrker       = Dir->make<TH2D>("hitsCAPHRIvTrker", "Number of CAPHRI Hits v Number of Tracker Hits; nTrackerHits; nCAPHRIHits", 100, 0.0, 7500, 100, 0.0, 100);

    //-------------1D Residuals-------------
    Hist->resid_nTimeClusters  = Dir->make<TH1D>("resid_nTimeClusters", "Residual between nPOT and predicted nPOT using the nTimeClusters; nPOT-PredictednPOT", 100, -1e6, 1e6);
    Hist->resid_nCaloHits      = Dir->make<TH1D>("resid_nCaloHits", "Residual between nPOT and predicted nPOT using the nCaloHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);
    Hist->resid_caloE_reco       = Dir->make<TH1D>("resid_caloE_reco", "Normalized Residual:  nPOT - reconstucted  nPOT from Calorimeter Energy; ( nPOT-recoNPOT)/nPOT", 100, -0.600, 0.600); //eventually -.5,5 i hope
    Hist->resid_nCaphriHits    = Dir->make<TH1D>("resid_nCAPHRIHits", "Residual between nPOT and predicted nPOT using the nCAPHRIHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);
    Hist->resid_nProtonTCs     = Dir->make<TH1D>("resid_nProtonTCs", "Residual between nPOT and predicted nPOT using the nProtonTCs;  nPOT-PredictednPOT", 100, -1e8, 1e8);
    Hist->resid_nTrackerHits   = Dir->make<TH1D>("resid_nTrackerHits", "Residual between nPOT and predicted nPOT using the nTrackerHits;  nPOT-PredictednPOT", 100, -1e7, 1e7);


    //------------Residual v nPOT-------------
    Hist->resid_nTimeClusters2d = Dir->make<TH2D>("resid_nTimeClusters2d", "nPOT v Residual of nPOT from nTimeClusters; nPOT-PredictednPOT; nPOT", 100, -1e6, 1e6, 100, 0.0, 1e8);
    Hist->resid_nCaloHits2d     = Dir->make<TH2D>("resid_nCaloHits2d", "nPOT v Residual of nPOT from nCaloHits; nPOT-PredictednPOT; ", 100, -1e7, 1e7, 100, 0.0, 1e8);
    Hist->resid_caloE_reco2d       = Dir->make<TH2D>("resid_caloE_reco2d", "Normalized Residual of nPOT from Calorimeter Energy v nPOT;nPOT; (nPOT-recoNPOT)/nPOT", 100, 0.0, 1e8, 100, -0.600, 0.600);
    Hist->resid_nCaphriHits2d    = Dir->make<TH2D>("resid_nCAPHRIHits2d", "nPT v Residual of nPOT from nCAPHRIHits; nPOT-PredictednPOT; nPOT", 100, -1e7, 1e7, 100, 0.0, 1e8);
    Hist->resid_nProtonTCs2d    = Dir->make<TH2D>("resid_nProtonTCs2d", "nPOT v Residual of nPOT from nProtonTCs; nPOT-PredictednPOT; nPOT", 100, -1e8, 1e8, 100, 0.0, 1e8);
    Hist->resid_nTrackerHits2d  = Dir->make<TH2D>("resid_nTrackerHits2d", "nPOT v Residual of nPOT from nTrackerHits; nPOT-PredictednPOT; nPOT", 100, -1e7, 1e7, 100, 0.0, 1e8);

    //-------------Residual v Observable------------
    Hist->resid_caloE_reco_v_caloE2d = Dir->make<TH2D>("resid_caloE_reco_v_caloE2d", "Normalized Residual of nPOT from Calorimeter Energy v Calorimeter Energy; CaloE [MeV]; (nPOT-recoNPOT)/nPOT", 100, 0.0, 1e4, 100, -0.600, 0.600);
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

  }

  //-----------------------------------------------------------------------------
  void NPOTAnalysis::fillEventHistograms(EventHists* Hist) {

    //----------------1D Observables----------------
    Hist->nTimeClusters->Fill(_eventData.nClusters_);
    Hist->nCaloHits->Fill(_eventData.nCaloHits_);
    Hist->caloEnergy->Fill(_eventData.caloEnergy_);
    Hist->nCaphriHits->Fill(_eventData.nCaphriHits_);
    Hist->nProtonTCs->Fill(_eventData.nProtonTCs_);
    Hist->nTrackerHits->Fill(_eventData.nTrackerHits_);
    Hist->nProtonsDeltaFinder->Fill(_eventData.nProtonsDeltaFinder_);

    //----------------MC Truth nPOT----------------
    Hist->nPOT->Fill(_eventData.nPOT_);

    //----------------Reconstructed nPOT----------------
    Hist->recoNPOT_caloE->Fill(_eventData.recoNPOT_caloE_);

    //----------------Passed Filter----------------
    if (_eventData.passed_) {Hist->filteredRecoNPOT->Fill(_eventData.recoNPOT_caloE_); }

    //----------------nPOT v Observable----------------
    Hist->nTimeClusters2d->Fill(_eventData.nClusters_, _eventData.nPOT_);
    Hist->nCaloHits2d->Fill(_eventData.nCaloHits_, _eventData.nPOT_);
    Hist->caloEnergy2d->Fill(_eventData.caloEnergy_, _eventData.nPOT_);
    Hist->nCaphriHits2d->Fill(_eventData.nCaphriHits_, _eventData.nPOT_);
    Hist->nProtonTCs2d->Fill(_eventData.nProtonTCs_, _eventData.nPOT_);
    Hist->nTrackerHits2d->Fill(_eventData.nTrackerHits_, _eventData.nPOT_);
    Hist->nPOTvnProtonTracks->Fill(_eventData.nProtonsDeltaFinder_, _eventData.nPOT_);

    //---------------- nProtonsDeltaFinder v Observable----------------
    Hist->nTimeClusters2df->Fill(_eventData.nClusters_, _eventData.nProtonsDeltaFinder_);
    Hist->nCaloHits2df->Fill(_eventData.nCaloHits_, _eventData.nProtonsDeltaFinder_);
    Hist->caloEnergy2df->Fill(_eventData.caloEnergy_, _eventData.nProtonsDeltaFinder_);
    Hist->nCaphriHits2df->Fill(_eventData.nCaphriHits_, _eventData.nProtonsDeltaFinder_);
    Hist->nProtonTCs2df->Fill(_eventData.nProtonTCs_, _eventData.nProtonsDeltaFinder_);
    Hist->nTrackerHits2df->Fill(_eventData.nTrackerHits_, _eventData.nProtonsDeltaFinder_);

    //----------------Observable v Observable----------------
    Hist->caloEH->Fill(_eventData.nCaloHits_, _eventData.caloEnergy_);
    Hist->cTCs->Fill(_eventData.nProtonTCs_, _eventData.nClusters_);
    Hist->hitsTC->Fill(_eventData.nCaloHits_, _eventData.nTrackerHits_);
    Hist->hitsCAPHRIvCal->Fill(_eventData.nCaloHits_, _eventData.nCaphriHits_);
    Hist->hitsCAPHRIvTrker->Fill(_eventData.nTrackerHits_, _eventData.nCaphriHits_);


    //----------------1D Residual----------------
    Hist->resid_nTimeClusters->Fill(_eventData.resid_nTimeClusters_);
    Hist->resid_nCaloHits->Fill(_eventData.resid_nCaloHits_);
    Hist->resid_caloE_reco->Fill(_eventData.resid_caloE_reco_);
    Hist->resid_nCaphriHits->Fill(_eventData.resid_nCaphriHits_);
    Hist->resid_nProtonTCs->Fill(_eventData.resid_nProtonTCs_);
    Hist->resid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_);

    //----------------Residual v nPOT----------------
    Hist->resid_nTimeClusters->Fill(_eventData.resid_nTimeClusters_, _eventData.nPOT_);
    Hist->resid_nCaloHits->Fill(_eventData.resid_nCaloHits_, _eventData.nPOT_);
    Hist->resid_caloE_reco2d->Fill( _eventData.nPOT_, _eventData.resid_caloE_reco_); //residuals on y-axis
    Hist->resid_nCaphriHits->Fill(_eventData.resid_nCaphriHits_, _eventData.nPOT_);
    Hist->resid_nProtonTCs->Fill(_eventData.resid_nProtonTCs_, _eventData.nPOT_);
    Hist->resid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_, _eventData.nPOT_);

    //----------------Residual v Observable----------------
    Hist->resid_caloE_reco_v_caloE2d->Fill(_eventData.caloEnergy_, _eventData.resid_caloE_reco_);
  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::retrieveData() {

    // reset information
    _eventData.reset();

    //--------------Get if passed RecoNPOTFilter_module.cc--------------
    const art::InputTag triggerTag{"TriggerResults::" + _processName};
    art::Handle<art::TriggerResults> trigResultsH;
    _event->getByLabel(triggerTag, trigResultsH);
    if (trigResultsH.isValid()) {
      const art::TriggerResults* trigResults = trigResultsH.product();
      TriggerResultsNavigator trigNavig(trigResults);
      _eventData.passed_ = trigNavig.accepted(_filterName);
    } else {
      std::cout << "[NPOTAnalysis::" << __func__ << "] TriggerResults handle is not valid!" << std::endl;
      _eventData.passed_ = false;
    }

    //--------------Get Reconstructed nPOT from RecoNPOTMaker_module.cc--------------
    art::Handle<RecoProtonBunchIntensity> recoNPOTH;
    _event->getByLabel(_recoPOTTag,recoNPOTH);
    if(recoNPOTH.isValid()) {
      _eventData.recoNPOT_caloE_ = recoNPOTH->intensity();
    }
    else { std::cout << "[NPOTAnalysis::" << __func__ << "] Could not retrieve RecoNPOT!" << std::endl; }

    //--------------Get the EventData from ProtonBunchIntensity--------------
    art::Handle<ProtonBunchIntensity> MCPOTH;
    _event->getByLabel(_MCPOTTag, MCPOTH);
    if (MCPOTH.isValid()) {_eventData.nPOT_  = MCPOTH->intensity();}

    //--------------Get the EventData from IntensityInfoCalo--------------
    art::Handle<IntensityInfoCalo> caloInt;
    _event->getByLabel(_caloTag, caloInt);
    if (caloInt.isValid()){
      _eventData.nCaloHits_   = caloInt->nCaloHits();
      _eventData.caloEnergy_  = caloInt->caloEnergy();
      _eventData.nCaphriHits_ = caloInt->nCaphriHits();
    } else { std::cout << "Calorimeter intensity info is not valid" << std::endl;}

    //--------------Get the EventData from IntensityInfoTimeCluster--------------
    //Default is trigger-level TZClusterFinder
    art::Handle<IntensityInfoTimeCluster> tcInt;
    _event->getByLabel(_tcCollTag, tcInt);
    if (tcInt.isValid()) { _eventData.nProtonTCs_  = tcInt->nProtonTCs();}
    else { std::cout << "Time cluster intensity info is not valid" << std::endl;}

    //Use DeltaFinder Module
    art::Handle<IntensityInfoTimeCluster> dfInt;
    _event->getByLabel(_dfCollTag, dfInt);
    if (dfInt.isValid()) { _eventData.nProtonsDeltaFinder_ = dfInt->nProtonTCs();}
    else { std::cout << "Delta-Finder intensity info is not valid" << std::endl;}

    //--------------Get the EventData from IntensityInfoTrackerHits--------------
    art::Handle<IntensityInfoTrackerHits> thInt;
    _event->getByLabel(_thCollTag, thInt);
    if (thInt.isValid()) { _eventData.nTrackerHits_ = thInt->nTrackerHits();}
    else { std::cout << "Tracker hits intensity info is not valid" << std::endl;}

    //--------------Get the TimeClusterCollection--------------
    art::Handle<mu2e::TimeClusterCollection> tcCollHandle;
    _event->getByLabel(_tcCollTag, tcCollHandle);
    if (tcCollHandle.isValid()) {_tcColl = tcCollHandle.product();}
    else {_tcColl = nullptr;}
  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::computeEventData() {

    // compute event level information
    _eventData.nClusters_ = (_tcColl) ? (int)_tcColl->size() : -1;

    if (_eventData.nPOT_ > 0 && _eventData.recoNPOT_caloE_ > 0) {
      const double residual = _eventData.nPOT_ - _eventData.recoNPOT_caloE_;
      _eventData.resid_caloE_reco_ = residual / _eventData.nPOT_;
    } else {
      _eventData.resid_caloE_reco_ = 0.;
      std::cout << "Failed to calculate normalized residual: nPOT = "
                << _eventData.nPOT_ << ", recoNPOT = " << _eventData.recoNPOT_caloE_ << std::endl;
    }
  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::analyze(const art::Event& event) {

    // get the event
    _event = &event;

    // retrieve the event data
    retrieveData();

    // process the data
    computeEventData();

    // fill histograms
    fillHistograms();
  }
}
DEFINE_ART_MODULE(mu2e::NPOTAnalysis)
