//////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
// this module doesn't do reconstruction
// on input, it takes a list  of StrawHitFlags flags and evaluates performance
// of the delta electron tagging
//
// hit type = 0 : proton
//            1 : ele 0 < p < 20
//            2 : ele 20 < p < 90
//            3 : ele 90 < p < 110
//            4 : everything else
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TVector2.h"
// data
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
// Utilities
#include "Offline/Mu2eUtilities/inc/SimParticleTimeOffset.hh"
// diagnostics

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "Offline/CalPatRec/inc/HlPrint.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class DeltaFinderAna : public art::EDProducer {

  public:
    enum { kNStations      = 20 };
    enum { kNFaces         =  4 };
    enum { kNPanelsPerFace =  3 };
    enum { kc2_cut         = 50 };
    enum { kPerpRes        = 10 };
    enum { kMaxDxy         = 50 };
    enum { kmax_gap        = 1  };

    enum {
      kNEventHistSets    =  10,
      kNStrawHitHistSets = 100,
      kNMcHistSets       = 200
    };

  protected:

    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom;
      TH1F*  fNHits;
      TH2F*  fNHitsVsMom;
      TH2F*  fNHitsRVsMom;
      TH2F*  fNHitsRVsNHits;
      TH1F*  fNHitsDelta;
      TH1F*  fFractReco;
      TH2F*  fFractRecoVsNHits;
    };

    struct StrawHitHist_t {
      TH1F*  fType;
      TH1F*  fTime;
      TH1F*  fMom;                        // momentum of the particle which produced the hit
      TH1F*  fEnergyDep;
      TH1F*  fDeltaT;
      TH1F*  fPDGCode;
      TH1F*  fDeltaFlag;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNSecondHits;
      TH1F*  fNSecondHitsT;
      TH1F*  fNMc;
      TH1F*  fPDGCode;
      TH1F*  fNMcHits;
      TH1F*  fNHitsDeltaT;
      TH1F*  fNHitsDeltaR;
    };

    struct Hist_t {
      EventHist_t*    fEvent   [kNEventHistSets   ];
      StrawHitHist_t* fStrawHit[kNStrawHitHistSets];
      McHist_t*       fMc      [kNMcHistSets      ];
    };

    Hist_t  _hist;

//-----------------------------------------------------------------------------
// diagnostics structures
//-----------------------------------------------------------------------------
    struct McPart_t {
      int   fFirstStation;
      int   fLastStation;
      int   fNHitsDelta;            // number of hits associated with all reconstructed delta electrons
      int   fTime;                  // lowest out of the hit times

      const SimParticle*               fSim;
//-----------------------------------------------------------------------------
// a hit and its flag here could be inconsistent - the flag comes from the 'combo' combohit
//-----------------------------------------------------------------------------
      std::vector<const ComboHit*>     fListOfHits; // 1-straw combo hits
      std::vector<const StrawHitFlag*> fListOfFlags;// a parallel list, but the flags are those of "Combo" ComboHits

      McPart_t(const SimParticle* Sim = NULL) {
        fSim          = Sim;
        fFirstStation = 999;
        fLastStation  = -1;
        fNHitsDelta   = 0;
        fTime         = 1.e6;
      }

      ~McPart_t() {
      }

      const ComboHit*     Hit(int I)   const { return fListOfHits[I] ;    }
      int                 NHits()      const { return fListOfHits.size(); }
      int                 NHitsDelta() const { return fNHitsDelta;        }
      const StrawHitFlag* Flag (int I) const { return fListOfFlags[I];    }

      float Momentum() const {
        float px = fSim->startMomentum().px();
        float py = fSim->startMomentum().py();
        float pz = fSim->startMomentum().pz();
        return sqrt(px*px+py*py+pz*pz);
      }

      float Time() const { return fTime; }

    };

    struct McHitInfo_t {
      const  McPart_t*     fMc;
      const  StrawHitFlag* fFlag;
      int                  fType;   // 0:p, 1:ele p<20, 2:ele 20<p<80  3:ele 100<p<110 4:everything else
    };

//-----------------------------------------------------------------------------
// NStations stations, 4-1=3 faces (for hit w/ lower z), 3 panels (for hit w/ lower z)
// 2017-07-27 P.Murat: the 2nd dimension should be 3, right?
//-----------------------------------------------------------------------------
    std::vector<McPart_t*>    _list_of_mc_particles; // list_of_particles with hits in the tracker
    std::vector<McHitInfo_t>  _list_of_mc_hit_info ; // for each 1-straw hit, pointer to the MC info
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag                  _chCollTag;
    art::InputTag                  _shCollTag;              // straw hits        by "makeSH"
    art::InputTag                  _schCollTag;             // 1-straw combohits by "makeSH"
    art::InputTag                  _shfCollTag;
    art::InputTag                  _sdmcCollTag;
    int                            _debugLevel;
    int                            _diagLevel;
    int                            _printElectrons;         //
    int                            _printElectronsMinNHits;
    float                          _printElectronsMaxFReco;
    int                            _printElectronHits;
//-----------------------------------------------------------------------------
// cache of event or geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl;
    const ComboHitCollection*      _schColl; // one straw hit per combo hit
    const StrawHitCollection*      _shColl;
    const StrawHitFlagCollection*  _shfColl;
    const StrawDigiMCCollection*   _sdmcColl;

    const Tracker*                 _tracker;
    int                            _eventNum;
    int                            _nComboHits;
    int                            _nStrawHits; // not the total number of straw hits, but the
                                                // number of straw hits from combohits

    int                            fNHitsDeltaTot;
    int                            fNHitsDeltaReco;

    HlPrint*                       _hlp;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit     DeltaFinderAna(fhicl::ParameterSet const&);
    virtual      ~DeltaFinderAna();

    void bookEventHistograms   (EventHist_t*    Hist, int HistSet, art::TFileDirectory* Dir);
    void bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void bookMcHistograms      (McHist_t*       Hist, int HistSet, art::TFileDirectory* Dir);

    void bookHistograms();

    void fillEventHistograms   (EventHist_t*    Hist);
    void fillMcHistograms      (McHist_t*       Hist, McPart_t* Mc );
    void fillStrawHitHistograms(StrawHitHist_t* Hist, const StrawHit* Hit, McHitInfo_t* McHitInfo);

    void fillHistograms();

    //    void orderHits();

    // void printComboHit(const ComboHit* Hit, const StrawGasStep* Step,
    //                    const char* Opt, int IHit, int Flags);
    void debug();
//-----------------------------------------------------------------------------
// overloaded methods of the base class
//-----------------------------------------------------------------------------
    virtual void beginJob();
    virtual void beginRun(art::Run& ARun);
    virtual void produce( art::Event& e);

//-----------------------------------------------------------------------------
// do these need to be private ?
//-----------------------------------------------------------------------------
    bool      findData     (const art::Event&  Evt);
    McPart_t* findParticle (const SimParticle* Sim);
    int       initMcDiag      ();
    int       associateMcTruth();
  };

//-----------------------------------------------------------------------------
  DeltaFinderAna::DeltaFinderAna(fhicl::ParameterSet const& pset):
    art::EDProducer(pset),
    _chCollTag             (pset.get<string>       ("chCollTag"             )),
    _shCollTag             (pset.get<string>       ("shCollTag"             )),
    _schCollTag            (pset.get<string>       ("schCollTag"            )),
    _shfCollTag            (pset.get<string>       ("shfCollTag"            )),
    _sdmcCollTag           (pset.get<art::InputTag>("sdmcCollTag"           )),
    _debugLevel            (pset.get<int>          ("debugLevel"            )),
    _diagLevel             (pset.get<int>          ("diagLevel"             )),
    _printElectrons        (pset.get<int>          ("printElectrons"        )),
    _printElectronsMinNHits(pset.get<int>          ("printElectronsMinNHits")),
    _printElectronsMaxFReco(pset.get<float>        ("printElectronsMaxFReco")),
    _printElectronHits     (pset.get<int>          ("printElectronHits"     ))
  {
    _hlp = HlPrint::Instance();
  }

  DeltaFinderAna::~DeltaFinderAna() {
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookEventHistograms(EventHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>(Form("event_%02i", HistSet), "Event Number", 100, 0., 100000.);

    Hist->fNSecondHits     = Dir->make<TH1F>(Form("nhit2_%02i" ,HistSet), "N(second hits)", 100, 0., 100.);
    Hist->fNSecondHitsT    = Dir->make<TH1F>(Form("nhit2t_%02i",HistSet), "N(second hits) w/ time selection", 100, 0., 100.);

    Hist->fNMc             = Dir->make<TH1F>(Form("nmc_%02i"      ,HistSet), "N(MC particles)", 100, 0., 1000.);
    Hist->fPDGCode         = Dir->make<TH1F>(Form("pdg_code_%02i" ,HistSet), "PDG Code"       , 100, 0., 100.);
    Hist->fNMcHits         = Dir->make<TH1F>(Form("n_mc_hits_%02i",HistSet), "N(MC hits)"     , 250, 0., 500.);
    Hist->fNHitsDeltaT     = Dir->make<TH1F>(Form("n_delta_ht_%02i" ,HistSet), "N(delta hits T)", 500, 0., 5000.);
    Hist->fNHitsDeltaR     = Dir->make<TH1F>(Form("n_delta_hr_%02i" ,HistSet), "N(delta hits R)", 500, 0., 5000.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookMcHistograms(McHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fPDGCode    = Dir->make<TH1F>("pdg"         , "PDG code"        , 500, -250., 250.);
    Hist->fMom        = Dir->make<TH1F>("mom"         , "momentum"        , 200, 0., 200.);
    Hist->fNHits      = Dir->make<TH1F>("nhits"       , "N(hits)"         , 200, 0., 200.);
    Hist->fNHitsVsMom = Dir->make<TH2F>("nhits_vs_mom", "N(hits) vs mom"  , 100,0,100, 200, 0., 200.);
    Hist->fNHitsRVsMom = Dir->make<TH2F>("nhr_vs_mom", "N(hits R) vs mom"  , 100,0,100, 200, 0., 200.);
    Hist->fNHitsRVsNHits = Dir->make<TH2F>("nhr_vs_nh", "N(hits R) vs NH"  , 100,0,100, 100, 0., 100.);
    Hist->fNHitsDelta = Dir->make<TH1F>("nhitsr"      , "N(hits reco)"   , 200, 0., 200.);
    Hist->fFractReco  = Dir->make<TH1F>("fractr"      , "NR/N"           , 100, 0.,   1.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fType      = Dir->make<TH1F>("type"  , "Hit type"        ,   10, 0., 10.);
    Hist->fTime      = Dir->make<TH1F>("time"  , "time"            ,  400, 0., 2000.);
    Hist->fMom       = Dir->make<TH1F>("mom"   , "Momentum"        ,  800, 0., 400.);
    Hist->fEnergyDep = Dir->make<TH1F>("edep"  , "edep"            ,  200, 0., 2e-2);
    Hist->fDeltaT    = Dir->make<TH1F>("dt"    , "DeltaT"          ,  200, -10,10);
    Hist->fPDGCode   = Dir->make<TH1F>("pdg"   , "PDG code"        , 2000, -10000,10000);
    Hist->fDeltaFlag = Dir->make<TH1F>("deflag", "Delta Flag = 1"  ,    2, 0,1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookHistograms() {
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    art::ServiceHandle<art::TFileService> tfs;
    char   folder_name[200];

    TH1::AddDirectory(0);

    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],i,&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book straw hit histograms
//-----------------------------------------------------------------------------
    int book_straw_hit_histset[kNStrawHitHistSets];
    for (int i=0; i<kNStrawHitHistSets; i++) book_straw_hit_histset[i] = 0;

    book_straw_hit_histset[  0] = 1;                // all
    book_straw_hit_histset[  1] = 1;                // all prot and deut
    book_straw_hit_histset[  2] = 1;                // all e-: p<20
    book_straw_hit_histset[  3] = 1;                // all e- 20<p<80
    book_straw_hit_histset[  4] = 1;                // all e-: 80<p<110
    book_straw_hit_histset[  5] = 1;                // all e-  p > 110
    book_straw_hit_histset[  6] = 1;                // all e+
    book_straw_hit_histset[  7] = 1;                // all mu- and mu+
    book_straw_hit_histset[  8] = 1;                // all everything else

    book_straw_hit_histset[ 10] = 1;                // delta
    book_straw_hit_histset[ 11] = 1;                // delta prot and deut
    book_straw_hit_histset[ 12] = 1;                // delta e-: p<20
    book_straw_hit_histset[ 13] = 1;                // delta e- 20<p<80
    book_straw_hit_histset[ 14] = 1;                // delta e-: 80<p<110
    book_straw_hit_histset[ 15] = 1;                // delta e-  p > 110
    book_straw_hit_histset[ 16] = 1;                // delta e+
    book_straw_hit_histset[ 17] = 1;                // delta mu- and mu+
    book_straw_hit_histset[ 18] = 1;                // delta everything else

    book_straw_hit_histset[ 20] = 1;                // non-delta
    book_straw_hit_histset[ 21] = 1;                // non-delta prot and deut
    book_straw_hit_histset[ 22] = 1;                // non-delta e-: p<20
    book_straw_hit_histset[ 23] = 1;                // non-delta e- 20<p<80
    book_straw_hit_histset[ 24] = 1;                // non-delta e-: 80<p<110
    book_straw_hit_histset[ 25] = 1;                // non-delta e-  p > 110
    book_straw_hit_histset[ 26] = 1;                // non-delta e+
    book_straw_hit_histset[ 27] = 1;                // non-delta mu- and mu+
    book_straw_hit_histset[ 28] = 1;                // non-delta everything else

    for (int i=0; i<kNStrawHitHistSets; i++) {
      if (book_straw_hit_histset[i] != 0) {
        sprintf(folder_name,"sh_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fStrawHit[i] = new StrawHitHist_t;
        bookStrawHitHistograms(_hist.fStrawHit[i],i,&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book MC histograms
//-----------------------------------------------------------------------------
    int book_mc_histset[kNMcHistSets];
    for (int i=0; i<kNMcHistSets; i++) book_mc_histset[i] = 0;

    book_mc_histset[  0] = 1;                // all particles
    book_mc_histset[  1] = 1;                // electrons
    book_mc_histset[  2] = 1;                // electrons fTime > 550
    book_mc_histset[  3] = 1;                // electrons fTime > 550 with last>first
    book_mc_histset[  4] = 1;                // electrons fTime > 550 with last>first and 6+ hits
    book_mc_histset[  5] = 1;                // electrons fTime > 550 with last>first, 5+ hits, and reco delta
    book_mc_histset[  6] = 1;                // electrons fTime > 550 with last>first, 5+ hits and p < 20

    book_mc_histset[100] = 1;                // electrons fTime > 550 and 20 < p < 80 MeV/c
    book_mc_histset[101] = 1;                // electrons fTime > 550 and p > 80 MeV/c

    for (int i=0; i<kNMcHistSets; i++) {
      if (book_mc_histset[i] != 0) {
        sprintf(folder_name,"mc_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fMc[i] = new McHist_t;
        bookMcHistograms(_hist.fMc[i],i,&tfdir);
      }
    }
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::beginJob() {
    bookHistograms();
  }


//----Get data------------------------------------------------------------------------------------------------
  void DeltaFinderAna::beginRun(art::Run& aRun) {

    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    _tracker = ttHandle.get();
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillEventHistograms(EventHist_t* Hist) {
    Hist->fEventNumber->Fill(_eventNum);

//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      int n_mc_hits = mc->fListOfHits.size();

      Hist->fNMcHits->Fill(n_mc_hits);
      Hist->fPDGCode->Fill(mc->fSim->pdgId());
    }

    Hist->fNMc->Fill(nmc);
    Hist->fNHitsDeltaT->Fill(fNHitsDeltaTot);
    Hist->fNHitsDeltaR->Fill(fNHitsDeltaReco);
  }

//-----------------------------------------------------------------------------
// assume that, for a ComboHit,all hits are from the same particle
//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillStrawHitHistograms(StrawHitHist_t* Hist, const StrawHit* Hit, McHitInfo_t* McHitInfo) {

    const McPart_t* mc = McHitInfo->fMc;

    Hist->fType->Fill(McHitInfo->fType);
    Hist->fTime->Fill(Hit->time());
    Hist->fMom->Fill(mc->Momentum());
    Hist->fEnergyDep->Fill(Hit->energyDep());
    Hist->fDeltaT->Fill(Hit->dt());
    Hist->fPDGCode->Fill(mc->fSim->pdgId());
    Hist->fPDGCode->Fill(mc->fSim->pdgId());

    int delta_flag = McHitInfo->fFlag->hasAnyProperty(StrawHitFlag::bkg);
    Hist->fDeltaFlag->Fill(delta_flag);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::fillMcHistograms(McHist_t* Hist, McPart_t* Mc) {
    const SimParticle* sim = Mc->fSim;
    float mom = Mc->Momentum();

    Hist->fPDGCode->Fill(sim->pdgId());
    Hist->fMom->Fill(mom);
    Hist->fNHits->Fill(Mc->NHits());
    Hist->fNHitsVsMom->Fill(mom,Mc->NHits());
    Hist->fNHitsDelta->Fill(Mc->fNHitsDelta);
    Hist->fNHitsRVsMom->Fill(mom,Mc->NHitsDelta());
    Hist->fNHitsRVsNHits->Fill(Mc->NHits(),Mc->NHitsDelta());

    float freco = Mc->fNHitsDelta/(Mc->NHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNHits->Fill(Mc->NHits(),freco);
  }

//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillHistograms() {
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms, mc_hit_info relates to the straw hit type
// 0:p, 1:ele p<20, 2:ele 20<p<80  3:ele 100<p<110 4:everything else
//-----------------------------------------------------------------------------
    int loc = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit* ch         = &_chColl->at(i);
      int nsh = ch->nStrawHits();
      for (int ish=0; ish<nsh; ish++) {
        int ind = ch->indexArray().at(ish);

        const StrawHit* sh = &_shColl->at(ind);

        const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
        const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
        const SimParticle*  sim  = sgs->simParticle().get();

        float px = sim->startMomentum().x();
        float py = sim->startMomentum().y();
        float pz = sim->startMomentum().z();
        float momentum = sqrt(px*px+py*py+pz*pz);

        int mc_type = -1;

        if      ((sim->pdgId() == 2212) or (sim->pdgId() == 1000010020)) mc_type = 0;
        else if (sim->pdgId() == 11) {
          if      (momentum <  20) mc_type = 1;
          else if (momentum <  80) mc_type = 2;
          else if (momentum < 110) mc_type = 3;
          else                     mc_type = 4;
        }
        else if (sim->pdgId() == -11)                           mc_type = 5;
        else if ((sim->pdgId() == 13) or (sim->pdgId() == -13)) mc_type = 6;
        else                                                    mc_type = 7;

        McPart_t mc(sim);
        McHitInfo_t    mc_hit_info;

        mc_hit_info.fMc   = &mc;
        mc_hit_info.fType = mc_type;
        mc_hit_info.fFlag = &_shfColl->at(i);

        fillStrawHitHistograms(_hist.fStrawHit[0],sh,&mc_hit_info);  // all

        if      (mc_type == 0) fillStrawHitHistograms(_hist.fStrawHit[1],sh,&mc_hit_info);
        else if (mc_type == 1) fillStrawHitHistograms(_hist.fStrawHit[2],sh,&mc_hit_info);
        else if (mc_type == 2) fillStrawHitHistograms(_hist.fStrawHit[3],sh,&mc_hit_info);
        else if (mc_type == 3) fillStrawHitHistograms(_hist.fStrawHit[4],sh,&mc_hit_info);
        else if (mc_type == 4) fillStrawHitHistograms(_hist.fStrawHit[5],sh,&mc_hit_info);
        else if (mc_type == 5) fillStrawHitHistograms(_hist.fStrawHit[6],sh,&mc_hit_info);
        else if (mc_type == 6) fillStrawHitHistograms(_hist.fStrawHit[7],sh,&mc_hit_info);
        else                   fillStrawHitHistograms(_hist.fStrawHit[8],sh,&mc_hit_info);

        const StrawHitFlag* flag = mc_hit_info.fFlag;

        // int edepOK = flag->hasAllProperties(StrawHitFlag::energysel);
        int delta  = flag->hasAllProperties(StrawHitFlag::bkg);

        if (delta) {
//-----------------------------------------------------------------------------
// set 2: all hits flagged as delta electrons
//-----------------------------------------------------------------------------
          fillStrawHitHistograms(_hist.fStrawHit[10],sh,&mc_hit_info);

          if      (mc_type == 0) fillStrawHitHistograms(_hist.fStrawHit[11],sh,&mc_hit_info);
          else if (mc_type == 1) fillStrawHitHistograms(_hist.fStrawHit[12],sh,&mc_hit_info);
          else if (mc_type == 2) fillStrawHitHistograms(_hist.fStrawHit[13],sh,&mc_hit_info);
          else if (mc_type == 3) fillStrawHitHistograms(_hist.fStrawHit[14],sh,&mc_hit_info);
          else if (mc_type == 4) fillStrawHitHistograms(_hist.fStrawHit[15],sh,&mc_hit_info);
          else if (mc_type == 5) fillStrawHitHistograms(_hist.fStrawHit[16],sh,&mc_hit_info);
          else if (mc_type == 6) fillStrawHitHistograms(_hist.fStrawHit[17],sh,&mc_hit_info);
          else                   fillStrawHitHistograms(_hist.fStrawHit[18],sh,&mc_hit_info);
        }
        else {
          fillStrawHitHistograms(_hist.fStrawHit[20],sh,&mc_hit_info);

          if      (mc_type == 0) fillStrawHitHistograms(_hist.fStrawHit[21],sh,&mc_hit_info);
          else if (mc_type == 1) fillStrawHitHistograms(_hist.fStrawHit[22],sh,&mc_hit_info);
          else if (mc_type == 2) fillStrawHitHistograms(_hist.fStrawHit[23],sh,&mc_hit_info);
          else if (mc_type == 3) fillStrawHitHistograms(_hist.fStrawHit[24],sh,&mc_hit_info);
          else if (mc_type == 4) fillStrawHitHistograms(_hist.fStrawHit[25],sh,&mc_hit_info);
          else if (mc_type == 5) fillStrawHitHistograms(_hist.fStrawHit[26],sh,&mc_hit_info);
          else if (mc_type == 6) fillStrawHitHistograms(_hist.fStrawHit[27],sh,&mc_hit_info);
          else                   fillStrawHitHistograms(_hist.fStrawHit[28],sh,&mc_hit_info);
        }

        loc++;
      }
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      const SimParticle* sim = mc->fSim;

      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if (sim->pdgId() == PDGCode::e_minus) {
        fillMcHistograms(_hist.fMc[1],mc);
        if (mc->Time() > 550) {
          fillMcHistograms(_hist.fMc[2],mc);
          if (mc->fLastStation > mc->fFirstStation) {
            fillMcHistograms(_hist.fMc[3],mc);
            if (mc->NHits() >= 5) {
              fillMcHistograms(_hist.fMc[4],mc);

              if (mc->Momentum() < 20) fillMcHistograms(_hist.fMc[6],mc);
            }
//-----------------------------------------------------------------------------
// a closer look at misreconstructed delta electrons
//-----------------------------------------------------------------------------
            float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

            if ((mc->Momentum() < 5) && (mc->Time() > 550) && (mc->NHits() > 40) && (fr < 0.5)) {
              printf(" event: %6i missed delta: sim.id = %10li mom = %10.3f time= %9.3f nhits = %3i nhits(delta): %3i first: %2i last: %2i",
                     _eventNum,
                     sim->id().asInt(), mc->Momentum(), mc->Time(),
                     mc->NHits(), mc->fNHitsDelta,
                     mc->fFirstStation, mc->fLastStation);
              printf(" fraction: %6.3f\n",fr);
            }
          }

          if ((mc->Momentum() > 20) && (mc->Momentum() < 80)) fillMcHistograms(_hist.fMc[100],mc);
          if ( mc->Momentum() > 80)                           fillMcHistograms(_hist.fMc[101],mc);
        }
      }
    }
  }

//-----------------------------------------------------------------------------
  DeltaFinderAna::McPart_t* DeltaFinderAna::findParticle(const SimParticle* Sim) {
    McPart_t* found(0);

    int n = _list_of_mc_particles.size();

    for (int i=0; i<n; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      if (mc->fSim == Sim) {
        found = mc;
        break;
      }
    }

    return found;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker - hope, it is shorter
// than the list of all particles
//-----------------------------------------------------------------------------
  int DeltaFinderAna::initMcDiag() {

    int n = _list_of_mc_particles.size();
    for (int i=0; i<n; i++) {
      McPart_t* p = _list_of_mc_particles.at(i);
      delete p;
    }

    _list_of_mc_particles.clear();

    _list_of_mc_hit_info.clear();
//-----------------------------------------------------------------------------
// count the number of participating straw hits (remember, it could be less than
// the total number of straw hits :participating are only the straw hits from
// combo hits
//-----------------------------------------------------------------------------
    _nStrawHits = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch   = &_chColl->at(i);
      _nStrawHits += ch->nStrawHits();
    }

    _list_of_mc_hit_info.resize(_nStrawHits);

    int delta_nhits_tot = 0;
//-----------------------------------------------------------------------------
// 'shf' is the combohit flag, despite the name
// so some straw hits may have flags not consistent with their origin
//-----------------------------------------------------------------------------
    int loc = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch   = &_chColl->at(i);
      const StrawHitFlag* shf  = &_shfColl->at(i);

      int nsh = ch->nStrawHits();
      for (int ish=0; ish<nsh; ish++) {
                                        // 'ind' is not an index in _list_of_mc_hit_info ,
                                        // how to find it ? - count !

        int ind = ch->indexArray().at(ish);
        const ComboHit* sch = &_schColl->at(ind);
//-----------------------------------------------------------------------------
// MC truth is based on the first straw hit of the combohit
//-----------------------------------------------------------------------------
        const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
        const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
        const SimParticle*  sim  = sgs->simParticle().get();
//-----------------------------------------------------------------------------
// search if this particle has already been registered
//-----------------------------------------------------------------------------
        McPart_t* mc = findParticle(sim);

        if (mc == NULL) {
                                        // add new particle
          mc = new McPart_t(sim);
          _list_of_mc_particles.push_back(mc);
        }
        mc->fListOfHits.push_back(sch);

        int station = sch->strawId().getStation();

        if (station < mc->fFirstStation) mc->fFirstStation = station;
        if (station > mc->fLastStation ) mc->fLastStation  = station;

        if (sch->correctedTime() < mc->fTime) mc->fTime = sch->correctedTime();
//-----------------------------------------------------------------------------
// remember, this is a combohit flag !
//-----------------------------------------------------------------------------
        mc->fListOfFlags.push_back(shf);

        McHitInfo_t* mc_hit_info = &_list_of_mc_hit_info.at(loc);

        mc_hit_info->fMc   = mc;
        mc_hit_info->fFlag = shf;

        int pdg_id = mc->fSim->pdgId();

        if      (pdg_id == 2212)   mc_hit_info->fType = 0;
        else if (pdg_id == 11  ) {
          float mom = mc->Momentum();
          if      (mom <  20)      {
            mc_hit_info->fType = 1;
            delta_nhits_tot++;
          }
          else if (mom <  90)       mc_hit_info->fType = 2;
          else if (mom < 110)       mc_hit_info->fType = 3;
          else                      mc_hit_info->fType = 4;
        }
        else if (pdg_id      == -11)mc_hit_info->fType = 5;
        else if (abs(pdg_id) == 13) mc_hit_info->fType = 6;
        else                        mc_hit_info->fType = 7;

        int flagged_as_delta = shf->hasAnyProperty(StrawHitFlag::bkg);

        if (flagged_as_delta) fNHitsDeltaReco++;
        loc++;
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// form a list of MC particles with hits in the tracker
//-----------------------------------------------------------------------------
  int DeltaFinderAna::associateMcTruth() {
//-----------------------------------------------------------------------------
// for each MC electron calculate the number of reconstructed hits
//-----------------------------------------------------------------------------
    int nmc    = _list_of_mc_particles.size();

    StrawHitFlag        deltamask(StrawHitFlag::bkg);

    for (int i=0; i<nmc; i++) {
      McPart_t* mc    = _list_of_mc_particles.at(i);
      mc->fNHitsDelta = 0;
//-----------------------------------------------------------------------------
// loop over the hits of MC delta electron and calculate fraction of them
// which have been tagged as the delta electron hits
//-----------------------------------------------------------------------------
      int nh = mc->fListOfHits.size();
      for (int ih=0; ih<nh; ih++) {
        const StrawHitFlag* flag = mc->fListOfFlags.at(ih);

        int flagged_as_delta = flag->hasAnyProperty(deltamask);

        if (flagged_as_delta) mc->fNHitsDelta += 1;
      }

      int pdg_id = mc->fSim->pdgId();

      if (pdg_id == PDGCode::e_minus  ) {
        float mom = mc->Momentum();
        if (mom < 20)      {
//-----------------------------------------------------------------------------
// call this "a delta electron"
//-----------------------------------------------------------------------------
          fNHitsDeltaTot += mc->NHits();
        }
      }
    }

    //    printf("DeltaFinderAna::associateMcTruth: fNHitsDeltaTot = %5i\n",fNHitsDeltaTot);

    return 0;
  }

//-----------------------------------------------------------------------------
bool DeltaFinderAna::findData(const art::Event& Evt) {
    _chColl     = nullptr;
    _sdmcColl   = nullptr;
    _shfColl    = nullptr;

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();
    _nComboHits = _chColl->size();

    auto shcH   = Evt.getValidHandle<StrawHitCollection>(_shCollTag);
    _shColl     = shcH.product();
    //    _nStrawHits = _shColl->size();

    auto schcH  = Evt.getValidHandle<ComboHitCollection>(_schCollTag);
    _schColl    = schcH.product();

    auto shfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_shfCollTag);
    _shfColl    = shfcH.product();

    auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
    _sdmcColl   = sdmccH.product();

    return (_chColl != 0) && (_nComboHits > 0) && (_shfColl != 0) && (_sdmcColl != 0) ;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::produce(art::Event& Event) {

    _eventNum = Event.event();
    if (_debugLevel) {
      printf("* >>> DeltaFinderAna::%s event number: %10i\n",__func__,_eventNum);
    }

    _hlp->SetEvent(&Event);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      throw cet::exception("RECO")
        << "mu2e::DeltaFinderAna_module::produce: missing data" << endl;
    }

    fNHitsDeltaTot  = 0;
    fNHitsDeltaReco = 0;

    initMcDiag      ();
    associateMcTruth();
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    fillHistograms  ();

    if (_debugLevel    > 0) debug();
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  void DeltaFinderAna::debug() {
    // hlp->printComboHitCollection(_data->chCollTag.encode().data(),
    //                              _data->chfCollTag.encode().data(),
    //                              _data->sdmcCollTag.encode().data());
//-----------------------------------------------------------------------------
// print MC electrons
//-----------------------------------------------------------------------------
    if (_printElectrons) {
      int nmc = _list_of_mc_particles.size();

      for (int i=0; i<nmc; i++) {
        McPart_t* mc = _list_of_mc_particles.at(i);
        const SimParticle* sim = mc->fSim;

        if ((sim->pdgId() == PDGCode::e_minus) && (mc->NHits()  >= _printElectronsMinNHits)) {

          float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

          if (fr < _printElectronsMaxFReco) {

            printf("* electron: sim.id: %6li mom: %8.3f time: %8.2f nhits: %3i nhits(delta): %3i first: %2i last: %2i",
                   sim->id().asInt(), mc->Momentum(), mc->Time(),
                   mc->NHits(),
                   mc->fNHitsDelta,
                   mc->fFirstStation, mc->fLastStation);
            printf(" fraction: %6.3f\n",fr);
//-----------------------------------------------------------------------------
// check if want to print hits
//-----------------------------------------------------------------------------
            if (_printElectronHits != 0) {
                                        // these are the 1-straw ComboHits
              int nh  = mc->NHits();
              const ComboHit* sch0 = &(*_schColl)[0];
              _hlp->printComboHit(0,0,"banner",0,0);

              for (int i=0; i<nh; i++) {
                const ComboHit* sch = mc->Hit(i);
                                        // flags here are the combo hit flags, not 1-straw hit flags
                int flags  = *((int*) mc->Flag(i));
                int ind    = sch-sch0;

                const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
                const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
                _hlp->printComboHit(sch,sgs,"data",ind,flags);
              }
            }
          }
        }
      }
    }
  }

// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(DeltaFinderAna)

}
