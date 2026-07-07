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
      TH1D* nCaloHits;
      TH1D* nCaloHitsD0; //disk 0 calo hits
      TH1D* nCaloHitsD1; //disk 1 calo hits
      TH1D* caloEnergy;
      TH1D* nCaphriHits;
      TH1D* nTrackerHits;
      TH1D* nProtonsTZ; // proton time clusters from TZClusterFinder
      TH1D* nProtonsDF; // proton time clusters from DeltaFinder

      TH1D* nPOT; //MC Truth  nPOT

      TH1D* recoNPOT_caloE; //Produced in RecoNPOTMaker_module.cc in CalPatRec

      TH1D* filteredRecoNPOT; //Passed RecoNPOTFilter_module.cc

      //-------------Observable v. NPOT-------------
      TH2D* nPOTvnCaloHits;
      TH2D* nPOTvnCaloHitsD0;
      TH2D* nPOTvnCaloHitsD1;
      TH2D* nPOTvcaloEnergy;
      TH2D* nPOTvnCaphriHits;
      TH2D* nPOTvnTrackerHits;
      TH2D* nPOTvnProtonsTZ;
      TH2D* nPOTvnProtonsDF;

      //-------------Observable v Observable-------------
      TH2D* caloEH; //caloEnergy v nCaloHits
      TH2D* cTCs; //nTimeClusters v nProtonTCs (from IntensityCalo Time Clusters)
      TH2D* hitsTC; //nTrackerHits v nCaloHits
      TH2D* hitsCAPHRIvCal; //nCaphriHits v nCalHits
      TH2D* hitsCAPHRIvTrker; //nCaphriHits v nTrackerHits

      //--------------Residual histograms: (nPOT predicted from Observable-nPOT)/(nPOT)--------------
      //naming scheme: resid_Observable
      TH1D* resid_nCaloHits;
      TH1D* resid_nCaloHitsD0;
      TH1D* resid_nCaloHitsD1;
      TH1D* resid_caloEnergy;
      TH1D* resid_nCaphriHits;
      TH1D* resid_nTrackerHits;
      TH1D* resid_nProtonsTZ;
      TH1D* resid_nProtonsDF;

      //--------------nPOT v residual histograms:--------------
      TH2D* nPOTvresid_nCaloHits;
      TH2D* nPOTvresid_nCaloHitsD0;
      TH2D* nPOTvresid_nCaloHitsD1;
      TH2D* nPOTvresid_caloEnergy;
      TH2D* nPOTvresid_nCaphriHits;
      TH2D* nPOTvresid_nTrackerHits;
      TH2D* nPOTvresid_nProtonsTZ;
      TH2D* nPOTvresid_nProtonsDF;

      //-------------Observable v Residual-------------
      TH2D* resid_caloEnergy_v_caloEnergy;

    };

    struct Hists {
      EventHists* _EventHists[kNEventHistsSets];
    };

    struct eventData {
      long long nPOT_;
      //-------------Observables-------------
      int nCaloHits_;
      int nCaloHitsD0_;
      int nCaloHitsD1_;
      int caloEnergy_;
      int nCaphriHits_;
      int nTrackerHits_;
      int nProtonsTZ_;
      int nProtonsDF_;

      //-------------Made in RecoNPOTMaker_module.cc as <RecoProtonBunchIntensity>
      long long recoNPOT_caloE_;

      //-------------Residuals-------------
      double resid_nCaloHits_;
      double resid_nCaloHitsD0_;
      double resid_nCaloHitsD1_;
      double resid_caloEnergy_;
      double resid_nCaphriHits_;
      double resid_nTrackerHits_;
      double resid_nProtonsTZ_;
      double resid_nProtonsDF_;

      //-------------Passed RecoNPOTFilter_module.cc-------------
      bool passed_;

      eventData() { reset(); }

      void reset() {
        nPOT_ = 0;
        nCaloHits_ = 0;
        nCaloHitsD0_ = 0;
        nCaloHitsD1_ = 0;
        caloEnergy_ = 0;
        nCaphriHits_ = 0;
        nTrackerHits_ = 0;
        nProtonsTZ_ = 0;
        nProtonsDF_ = 0;
        recoNPOT_caloE_ = 0;
        resid_nCaloHits_ = 0;
        resid_nCaloHitsD0_ = 0;
        resid_nCaloHitsD1_ = 0;
        resid_caloEnergy_ = 0;
        resid_nCaphriHits_ = 0;
        resid_nTrackerHits_ = 0;
        resid_nProtonsTZ_ = 0;
        resid_nProtonsDF_ = 0;
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

    constexpr double max_pot = 2.e8;

    //-------------1D Observables-------------
    Hist->nCaloHits              = Dir->make<TH1D>("nCaloHits", "Number of Calorimeter Hits; N(calo hits)"                             , 100, 0.0, 5000);
    Hist->nCaloHitsD0            = Dir->make<TH1D>("nCaloHitsD0", "Number of Calorimeter Hits in Disk 0; N(hits)"                      , 100, 0.0, 3000);
    Hist->nCaloHitsD1            = Dir->make<TH1D>("nCaloHitsD1", "Number of Calorimeter Hits in Disk 1; N(hits)"                      , 100, 0.0, 2000);
    Hist->caloEnergy             = Dir->make<TH1D>("caloEnergy", "Calorimeter Energy; E[MeV]"                                          , 100, 0.0,12000);
    Hist->nCaphriHits            = Dir->make<TH1D>("nCaphriHits", "Number of CAPHRI Hits; nCAPHRIHits"                                 ,  30, 0.0,   30);
    Hist->nTrackerHits           = Dir->make<TH1D>("nTrackerHits", "Number of Tracker Hits from Intensity Information; N(tracker hits)", 100, 0.0, 8000);
    Hist->nProtonsTZ             = Dir->make<TH1D>("nProtonsTZ", "Number of Proton Time Clusters from TZClusterFinder; N(proton TCs)"  ,  40, 0.0,   40);
    Hist->nProtonsDF             = Dir->make<TH1D>("nProtonsDF", "Number of Proton Time Clusters from DeltaFinder; N(proton TCs)"      ,  40, 0.0,   40);

    //-------------MC Truth nPOT-------------
    Hist->nPOT                   = Dir->make<TH1D>("nPOT", "Number of Protons on Target; nPOT"                                         , 100, 0.0, max_pot);

    //-------------Reconstructed nPOT from Observable-------------
    Hist->recoNPOT_caloE         = Dir->make<TH1D>("recoNPOT_caloE", "Reconstructed nPOT from Calorimeter Energy; RecoNPOT from CaloE" , 100, 0.0, max_pot);

    //-------------Filtered-------------
    Hist->filteredRecoNPOT       = Dir->make<TH1D>("filteredRecoNPOT", "Reconstructed nPOT that passed Filter; Filtered ReconPOT"      , 100, 0.0, max_pot);

    //-------------Observable v nPOT-------------
    Hist->nPOTvnCaloHits         = Dir->make<TH2D>("nPOTvnCaloHits"   , "N(POT) vs N(Calo hits); N(calo hits); N(POT)"                 , 100, 0.0, 5000, 100, 0.0, max_pot);
    Hist->nPOTvnCaloHitsD0       = Dir->make<TH2D>("nPOTvnCaloHitsD0" , "N(POT) vs N(Calo hits) in disk 0; N(disk 0 calo hits); N(POT)", 100, 0.0, 4000, 100, 0.0, max_pot);
    Hist->nPOTvnCaloHitsD1       = Dir->make<TH2D>("nPOTvnCaloHitsD1" , "N(POT) vs N(Calo hits) in disk 1; N(disk 1 calo hits); N(POT)", 100, 0.0, 2000, 100, 0.0, max_pot);
    Hist->nPOTvcaloEnergy        = Dir->make<TH2D>("nPOTvcaloEnergy"  , "N(POT) vs Calorimeter Energy; E [MeV]; N(POT)"                , 100, 0.0,12000, 100, 0.0, max_pot);
    Hist->nPOTvnCaphriHits       = Dir->make<TH2D>("nPOTvnCaphriHits" , "N(POT) vs N(CAPHRI hits); N(CAPHRI hits); N(POT)"             ,  30, 0.0,   30, 100, 0.0, max_pot);
    Hist->nPOTvnTrackerHits      = Dir->make<TH2D>("nPOTvnTrackerHits", "N(POT) vs N(Tracker hits); N(tracker hits); N(POT)"           , 100, 0.0, 6000, 100, 0.0, max_pot);
    Hist->nPOTvnProtonsTZ        = Dir->make<TH2D>("nPOTvnProtonsTZ"  , "N(POT) vs N(TZ protons); N(protons); nPOT"                    ,  40, 0.0,   40, 100, 0.0, max_pot);
    Hist->nPOTvnProtonsDF        = Dir->make<TH2D>("nPOTvnProtonsDF"  , "N(POT) vs N(DF protons); N(protons); nPOT"                    ,  40, 0.0,   40, 100, 0.0, max_pot);

    //-------------Observable v Observable-------------

    Hist->caloEH                 = Dir->make<TH2D>("caloEH", "Calorimeter Energy vs the number of Calorimeter Hits; nCaloHits; E[MeV]", 100, 0.0, 1300, 100, 0.0, 7500);
    Hist->cTCs                   = Dir->make<TH2D>("cTCs", "Number of Time Clusters vs Number of Time Clusters from Intensity Information; nProtonTCs from IntensityInfo; nTCs", 100, 0.0, 100, 100, 0.0, 100);
    Hist->hitsTC                 = Dir->make<TH2D>("hitsTC", "Number of Tracker Hits v number of Calorimeter Hits; nCaloHits; nTrackerHits", 100.0, 0.0, 1300, 20, 0.0, 7500);
    Hist->hitsCAPHRIvCal         = Dir->make<TH2D>("hitsCAPHRIvCal", "Number of CAPHRI Hits v Number of Calorimeter Hits; nCaloHits; nCAPHRIHits", 100, 0.0, 1300, 100, 0.0, 100);
    Hist->hitsCAPHRIvTrker       = Dir->make<TH2D>("hitsCAPHRIvTrker", "Number of CAPHRI Hits v Number of Tracker Hits; nTrackerHits; nCAPHRIHits", 100, 0.0, 7500, 100, 0.0, 100);

    //-------------1D Residuals-------------
    Hist->resid_nCaloHits        = Dir->make<TH1D>("resid_nCaloHits"   , "Number of Calorimeter Hits residual;"                         , 100, -1., 1.);
    Hist->resid_nCaloHitsD0      = Dir->make<TH1D>("resid_nCaloHitsD0" , "Number of Calorimeter Hits in Disk 0 residual;"               , 100, -1., 1.);
    Hist->resid_nCaloHitsD1      = Dir->make<TH1D>("resid_nCaloHitsD1" , "Number of Calorimeter Hits in Disk 1 residual;"               , 100, -1., 1.);
    Hist->resid_caloEnergy       = Dir->make<TH1D>("resid_caloEnergy"  , "Calorimeter Energy residual;"                                 , 100, -1., 1.);
    Hist->resid_nCaphriHits      = Dir->make<TH1D>("resid_nCaphriHits" , "Number of CAPHRI Hits residual;"                              , 100, -1., 1.);
    Hist->resid_nTrackerHits     = Dir->make<TH1D>("resid_nTrackerHits", "Number of Tracker Hits from Intensity Information residual;"  , 100, -1., 1.);
    Hist->resid_nProtonsTZ       = Dir->make<TH1D>("resid_nProtonsTZ"  , "Number of Proton Time Clusters from TZClusterFinder residual;", 100, -1., 1.);
    Hist->resid_nProtonsDF       = Dir->make<TH1D>("resid_nProtonsDF"  , "Number of Proton Time Clusters from DeltaFinder residual;"    , 100, -1., 1.);


    //------------Residual v nPOT-------------
    Hist->nPOTvresid_nCaloHits   = Dir->make<TH2D>("nPOTvresid_nCaloHits"   , "N(POT) vs Number of Calorimeter Hits residual;"                         , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nCaloHitsD0 = Dir->make<TH2D>("nPOTvresid_nCaloHitsD0" , "N(POT) vs Number of Calorimeter Hits in Disk 0 residual;"               , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nCaloHitsD1 = Dir->make<TH2D>("nPOTvresid_nCaloHitsD1" , "N(POT) vs Number of Calorimeter Hits in Disk 1 residual;"               , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_caloEnergy  = Dir->make<TH2D>("nPOTvresid_caloEnergy"  , "N(POT) vs Calorimeter Energy residual;"                                 , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nCaphriHits = Dir->make<TH2D>("nPOTvresid_nCaphriHits" , "N(POT) vs Number of CAPHRI Hits residual;"                              , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nTrackerHits= Dir->make<TH2D>("nPOTvresid_nTrackerHits", "N(POT) vs Number of Tracker Hits from Intensity Information residual;"  , 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nProtonsTZ  = Dir->make<TH2D>("nPOTvresid_nProtonsTZ"  , "N(POT) vs Number of Proton Time Clusters from TZClusterFinder residual;", 100, -1., 1., 100, 0., 1e8);
    Hist->nPOTvresid_nProtonsDF  = Dir->make<TH2D>("nPOTvresid_nProtonsDF"  , "N(POT) vs Number of Proton Time Clusters from DeltaFinder residual;"    , 100, -1., 1., 100, 0., 1e8);

    //-------------Residual v Observable------------
    Hist->resid_caloEnergy_v_caloEnergy = Dir->make<TH2D>("resid_caloEnergy_v_caloEnergy",
                                                          "N(POT) residual from Calorimeter Energy v Calorimeter Energy; Energy [MeV]; (N(POT)- reco N(POT))/N(POT)",
                                                          100, 0.0, 1e4, 100, -0.600, 0.600);
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
    Hist->nCaloHits   ->Fill(_eventData.nCaloHits_   );
    Hist->nCaloHitsD0 ->Fill(_eventData.nCaloHitsD0_ );
    Hist->nCaloHitsD1 ->Fill(_eventData.nCaloHitsD1_ );
    Hist->caloEnergy  ->Fill(_eventData.caloEnergy_  );
    Hist->nCaphriHits ->Fill(_eventData.nCaphriHits_ );
    Hist->nTrackerHits->Fill(_eventData.nTrackerHits_);
    Hist->nProtonsTZ  ->Fill(_eventData.nProtonsTZ_  );
    Hist->nProtonsDF  ->Fill(_eventData.nProtonsDF_  );

    //----------------MC Truth nPOT----------------
    Hist->nPOT->Fill(_eventData.nPOT_);

    //----------------Reconstructed nPOT----------------
    Hist->recoNPOT_caloE->Fill(_eventData.recoNPOT_caloE_);

    //----------------Passed Filter----------------
    if (_eventData.passed_) {Hist->filteredRecoNPOT->Fill(_eventData.recoNPOT_caloE_); }

    //----------------nPOT v Observable----------------
    Hist->nPOTvnCaloHits   ->Fill(_eventData.nCaloHits_   , _eventData.nPOT_);
    Hist->nPOTvnCaloHitsD0 ->Fill(_eventData.nCaloHitsD0_ , _eventData.nPOT_);
    Hist->nPOTvnCaloHitsD1 ->Fill(_eventData.nCaloHitsD1_ , _eventData.nPOT_);
    Hist->nPOTvcaloEnergy  ->Fill(_eventData.caloEnergy_  , _eventData.nPOT_);
    Hist->nPOTvnCaphriHits ->Fill(_eventData.nCaphriHits_ , _eventData.nPOT_);
    Hist->nPOTvnTrackerHits->Fill(_eventData.nTrackerHits_, _eventData.nPOT_);
    Hist->nPOTvnProtonsTZ  ->Fill(_eventData.nProtonsTZ_  , _eventData.nPOT_);
    Hist->nPOTvnProtonsDF  ->Fill(_eventData.nProtonsDF_  , _eventData.nPOT_);

    //----------------Observable v Observable----------------
    // Hist->caloEH->Fill(_eventData.nCaloHits_, _eventData.caloEnergy_);
    // Hist->cTCs->Fill(_eventData.nProtonTCs_, _eventData.nClusters_);
    // Hist->hitsTC->Fill(_eventData.nCaloHits_, _eventData.nTrackerHits_);
    // Hist->hitsCAPHRIvCal->Fill(_eventData.nCaloHits_, _eventData.nCaphriHits_);
    // Hist->hitsCAPHRIvTrker->Fill(_eventData.nTrackerHits_, _eventData.nCaphriHits_);


    //----------------1D Residual----------------
    Hist->resid_nCaloHits   ->Fill(_eventData.resid_nCaloHits_   );
    Hist->resid_nCaloHitsD0 ->Fill(_eventData.resid_nCaloHitsD0_ );
    Hist->resid_nCaloHitsD1 ->Fill(_eventData.resid_nCaloHitsD1_ );
    Hist->resid_caloEnergy  ->Fill(_eventData.resid_caloEnergy_  );
    Hist->resid_nCaphriHits ->Fill(_eventData.resid_nCaphriHits_ );
    Hist->resid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_);
    Hist->resid_nProtonsTZ  ->Fill(_eventData.resid_nProtonsTZ_  );
    Hist->resid_nProtonsDF  ->Fill(_eventData.resid_nProtonsDF_  );

    //----------------Residual v nPOT----------------
    Hist->nPOTvresid_nCaloHits   ->Fill(_eventData.resid_nCaloHits_   , _eventData.nPOT_);
    Hist->nPOTvresid_nCaloHitsD0 ->Fill(_eventData.resid_nCaloHitsD0_ , _eventData.nPOT_);
    Hist->nPOTvresid_nCaloHitsD1 ->Fill(_eventData.resid_nCaloHitsD1_ , _eventData.nPOT_);
    Hist->nPOTvresid_caloEnergy  ->Fill(_eventData.resid_caloEnergy_  , _eventData.nPOT_);
    Hist->nPOTvresid_nCaphriHits ->Fill(_eventData.resid_nCaphriHits_ , _eventData.nPOT_);
    Hist->nPOTvresid_nTrackerHits->Fill(_eventData.resid_nTrackerHits_, _eventData.nPOT_);
    Hist->nPOTvresid_nProtonsTZ  ->Fill(_eventData.resid_nProtonsTZ_  , _eventData.nPOT_);
    Hist->nPOTvresid_nProtonsDF  ->Fill(_eventData.resid_nProtonsDF_  , _eventData.nPOT_);

    //----------------Residual v Observable----------------
    Hist->resid_caloEnergy_v_caloEnergy->Fill(_eventData.caloEnergy_, _eventData.resid_caloEnergy_);
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
      _eventData.nCaloHitsD0_ = caloInt->nCaloHitsD0();
      _eventData.nCaloHitsD1_ = caloInt->nCaloHitsD1();
      _eventData.caloEnergy_  = caloInt->caloEnergy();
      _eventData.nCaphriHits_ = caloInt->nCaphriHits();
    } else { std::cout << "Calorimeter intensity info is not valid" << std::endl;}

    //--------------Get the EventData from IntensityInfoTimeCluster--------------
    //Default is trigger-level TZClusterFinder
    art::Handle<IntensityInfoTimeCluster> tcInt;
    _event->getByLabel(_tcCollTag, tcInt);
    if (tcInt.isValid()) { _eventData.nProtonsTZ_  = tcInt->nProtonTCs();}
    else { std::cout << "Time cluster intensity info is not valid" << std::endl;}

    //Use DeltaFinder Module
    art::Handle<IntensityInfoTimeCluster> dfInt;
    _event->getByLabel(_dfCollTag, dfInt);
    if (dfInt.isValid()) { _eventData.nProtonsDF_ = dfInt->nProtonTCs();}
    else { std::cout << "Delta-Finder intensity info is not valid" << std::endl;}

    //--------------Get the EventData from IntensityInfoTrackerHits--------------
    art::Handle<IntensityInfoTrackerHits> thInt;
    _event->getByLabel(_thCollTag, thInt);
    if (thInt.isValid()) { _eventData.nTrackerHits_ = thInt->nTrackerHits();}
    else { std::cout << "Tracker hits intensity info is not valid" << std::endl;}
  }

  //--------------------------------------------------------------------------------
  void NPOTAnalysis::computeEventData() {

    // compute event level information
    if (_eventData.nPOT_ > 0 && _eventData.recoNPOT_caloE_ > 0) {
      const double residual = _eventData.nPOT_ - _eventData.recoNPOT_caloE_;
      _eventData.resid_caloEnergy_ = residual / _eventData.nPOT_;
    } else {
      _eventData.resid_caloEnergy_ = 0.;
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
