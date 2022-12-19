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
#include "TH2F.h"
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
      std::vector<const ComboHit*>     fListOfHits;
      std::vector<const StrawHitFlag*> fListOfFlags;

      McPart_t(const SimParticle* Sim = NULL) {
        fSim          = Sim;
        fFirstStation = 999;
        fLastStation  = -1;
        fNHitsDelta   = 0;
        fTime         = 1.e6;
      }

      ~McPart_t() {
      }

      const ComboHit* Hit(int I)   const { return fListOfHits.at(I) ; }
      int             NHits()      const { return fListOfHits.size(); }
      int             NHitsDelta() const { return fNHitsDelta;        }

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
    std::vector<McHitInfo_t>  _list_of_mc_hit_info ; // for each straw hit, pointer to the MC info
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag                  _chCollTag;
    art::InputTag                  _shCollTag;              // straw hits "makeSH"
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
    const StrawHitCollection*      _shColl;
    const StrawHitFlagCollection*  _shfColl;
    const StrawDigiMCCollection*   _sdmcColl;

    const Tracker*                 _tracker;
    int                            _eventNum;
    int                            _nComboHits;
    int                            _nStrawHits;

    int                            fNHitsDeltaTot;
    int                            fNHitsDeltaReco;
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

    void printComboHit(const ComboHit* Hit, const StrawGasStep* Step,
                       const char* Opt, int IHit, int Flags);
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
    _shfCollTag            (pset.get<string>       ("shfCollTag"            )),
    _sdmcCollTag           (pset.get<art::InputTag>("sdmcCollTag"           )),
    _debugLevel            (pset.get<int>          ("debugLevel"            )),
    _diagLevel             (pset.get<int>          ("diagLevel"             )),
    _printElectrons        (pset.get<int>          ("printElectrons"        )),
    _printElectronsMinNHits(pset.get<int>          ("printElectronsMinNHits")),
    _printElectronsMaxFReco(pset.get<float>        ("printElectronsMaxFReco")),
    _printElectronHits     (pset.get<int>          ("printElectronHits"     ))
  {
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

    Hist->fPDGCode    = Dir->make<TH1F>("pdg"  , "PDG code"        , 500, -250., 250.);
    Hist->fMom        = Dir->make<TH1F>("mom"  , "momentum"        , 200, 0., 200.);
    Hist->fNHits      = Dir->make<TH1F>("nhits", "N(hits)"         , 200, 0., 200.);
    Hist->fNHitsDelta = Dir->make<TH1F>("nhitsr", "N(hits reco)"   , 200, 0., 200.);
    Hist->fFractReco  = Dir->make<TH1F>("fractr", "NR/N"           , 100, 0.,   1.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fType      = Dir->make<TH1F>("type", "Hit type"        , 10, 0., 10.);
    Hist->fTime      = Dir->make<TH1F>("time", "time"            , 400, 0., 2000.);
    Hist->fMom       = Dir->make<TH1F>("mom" , "Momentum"        , 200, 0., 400.);
    Hist->fEnergyDep = Dir->make<TH1F>("edep", "edep"            , 200, 0., 2e-2);
    Hist->fDeltaT    = Dir->make<TH1F>("dt"  , "DeltaT"          , 200, -10,10);
    Hist->fPDGCode   = Dir->make<TH1F>("pdg" , "PDG code"        , 2000, -10000,10000);
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

    book_straw_hit_histset[  0] = 1;                // all hits
    book_straw_hit_histset[  1] = 1;                // all hits t>500
    book_straw_hit_histset[  2] = 1;                // all hits t>500
    book_straw_hit_histset[  3] = 1;                // all hits t>500 edepOK
    book_straw_hit_histset[  4] = 1;                // all hits t>500 edepOK non-delta
    book_straw_hit_histset[  5] = 1;                // all hits t>500 edepOK non-delta MC=delta (type 1)
    book_straw_hit_histset[  6] = 1;                // all hits t>500 edepOK delta
    book_straw_hit_histset[  7] = 1;                // all hits t>500 edepOK delta MC=delta (type 1)

    book_straw_hit_histset[ 10] = 1;                // all hits type=0
    book_straw_hit_histset[ 20] = 1;                // all hits type=1
    book_straw_hit_histset[ 30] = 1;                // all hits type=2
    book_straw_hit_histset[ 40] = 1;                // all hits type=3
    book_straw_hit_histset[ 50] = 1;                // all hits type=3

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
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::fillMcHistograms(McHist_t* Hist, McPart_t* Mc) {
    const SimParticle* sim = Mc->fSim;
    float mom = Mc->Momentum();

    Hist->fPDGCode->Fill(sim->pdgId());
    Hist->fMom->Fill(mom);
    Hist->fNHits->Fill(Mc->NHits());
    Hist->fNHitsDelta->Fill(Mc->fNHitsDelta);

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
// straw hit histograms, mc_hit_info relates to the straw hit
//-----------------------------------------------------------------------------
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit* ch          = &_chColl->at(i);
      McHitInfo_t*    mc_hit_info = &_list_of_mc_hit_info.at(i);

      int nsh = ch->nStrawHits();
      for (int ish=0; ish<nsh; ish++) {
        int ind = ch->indexArray().at(ish);
        const StrawHit* sh = &_shColl->at(ind);

        fillStrawHitHistograms(_hist.fStrawHit[0],sh,mc_hit_info);
        if (sh->time() > 500) {
          fillStrawHitHistograms(_hist.fStrawHit[1],sh,mc_hit_info);

          const StrawHitFlag* flag = mc_hit_info->fFlag;

          fillStrawHitHistograms(_hist.fStrawHit[2],sh,mc_hit_info);
          int edepOK = flag->hasAllProperties(StrawHitFlag::energysel);
          if (edepOK) {
            fillStrawHitHistograms(_hist.fStrawHit[3],sh,mc_hit_info);
            int delta = flag->hasAllProperties(StrawHitFlag::bkg);
            if (! delta) {
//-----------------------------------------------------------------------------
// StrawHit SET 4: hits not marked as delta electron hits
//          SET 5: hits of low energy electrons not marked as delta electron hits
//-----------------------------------------------------------------------------
              fillStrawHitHistograms(_hist.fStrawHit[4],sh,mc_hit_info);
              if (mc_hit_info->fType == 1) { // low-energy electrons
                fillStrawHitHistograms(_hist.fStrawHit[5],sh,mc_hit_info);
              }
            }
            else {
//-----------------------------------------------------------------------------
// StrawHit SET 6: hits marked as delta electron hits
//          SET 7: hits of low energy electrons marked as such
//-----------------------------------------------------------------------------
              fillStrawHitHistograms(_hist.fStrawHit[6],sh,mc_hit_info);
              if (mc_hit_info->fType == 1) { // low-energy electrons
                fillStrawHitHistograms(_hist.fStrawHit[7],sh,mc_hit_info);
              }
            }
          }
        }

        if (mc_hit_info->fType == 0) fillStrawHitHistograms(_hist.fStrawHit[10],sh,mc_hit_info);
        if (mc_hit_info->fType == 1) fillStrawHitHistograms(_hist.fStrawHit[20],sh,mc_hit_info);
        if (mc_hit_info->fType == 2) fillStrawHitHistograms(_hist.fStrawHit[30],sh,mc_hit_info);
        if (mc_hit_info->fType == 3) fillStrawHitHistograms(_hist.fStrawHit[40],sh,mc_hit_info);
        if (mc_hit_info->fType == 4) fillStrawHitHistograms(_hist.fStrawHit[50],sh,mc_hit_info);
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
          if (mc->Momentum()  > 80)                           fillMcHistograms(_hist.fMc[101],mc);
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

    _list_of_mc_hit_info.resize(_nComboHits);

    int delta_nhits_tot = 0;

    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch   = &_chColl->at(i);
      const StrawHitFlag* shf  = &_shfColl->at(i);
//-----------------------------------------------------------------------------
// MC truth is based on the first straw hit of the combohit
//-----------------------------------------------------------------------------
      int ind = ch->indexArray().at(0);
      const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = &(*sgs->simParticle());
//-----------------------------------------------------------------------------
// search if this particle has already been registered
//-----------------------------------------------------------------------------
      McPart_t* mc = findParticle(sim);

      if (mc == NULL) {
                                        // add new particle
        mc = new McPart_t(sim);
        _list_of_mc_particles.push_back(mc);
      }
      mc->fListOfHits.push_back(ch);

      const StrawHit* sh = &_shColl->at(ind);
      StrawId      shid  = sh->strawId();
      const Straw& straw = _tracker->getStraw(shid);
      int station = straw.id().getStation();
      if (station < mc->fFirstStation) mc->fFirstStation = station;
      if (station > mc->fLastStation ) mc->fLastStation  = station;

      if (ch->correctedTime() < mc->fTime) mc->fTime = ch->correctedTime();

      mc->fListOfFlags.push_back(shf);

      McHitInfo_t* mc_hit_info = &_list_of_mc_hit_info.at(i);

      mc_hit_info->fMc   = mc;
      mc_hit_info->fFlag = shf;

      int pdg_id = mc->fSim->pdgId();

      if      (pdg_id == PDGCode::proton)   mc_hit_info->fType = 0;
      else if (pdg_id == PDGCode::e_minus ) {
        float mom = mc->Momentum();
        if      (mom <  20)      {
          mc_hit_info->fType = 1;
          delta_nhits_tot++;
        }
        else if (mom <  90)      mc_hit_info->fType = 2;
        else if (mom < 110)      mc_hit_info->fType = 3;
      }
      else                       mc_hit_info->fType = 4;

      int flagged_as_delta = shf->hasAnyProperty(StrawHitFlag::bkg);

      if (flagged_as_delta) fNHitsDeltaReco++;
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
    _nStrawHits = _shColl->size();

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
      printf(">>> %s::%s event number: %10i\n",typeid(*this).name(), __func__,_eventNum);
    }
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

            printf(" electron: sim.id = %10li mom = %10.3f time= %9.3f nhits = %3i nhits(delta): %3i first: %2i last: %2i",
                   sim->id().asInt(), mc->Momentum(), mc->Time(),
                   mc->NHits(),
                   mc->fNHitsDelta,
                   mc->fFirstStation, mc->fLastStation);
            printf(" fraction: %6.3f\n",fr);
//-----------------------------------------------------------------------------
// check if want to print hits
//-----------------------------------------------------------------------------
            if (_printElectronHits != 0) {
              int nh  = mc->NHits();
              const ComboHit* ch0 = &(*_chColl)[0];
              printComboHit(0,0,"banner",0,0);
              for (int i=0; i<nh; i++) {
                const ComboHit* ch       = mc->Hit(i);
                int loc    = ch-ch0;
                int flags  = *((int*) &_shfColl->at(loc));

                int ind                  = ch->indexArray().at(0);
                const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
                const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
                printComboHit(ch,sgs,"data",loc,flags);
              }
            }
          }
        }
      }
    }
  }

//-----------------------------------------------------------------------------
// stolen from Stntuple/print/TAnaDump.cc
//-----------------------------------------------------------------------------
void DeltaFinderAna::printComboHit(const ComboHit*     Hit,
                                   const StrawGasStep* Step,
                                   const char*         Opt, int IHit, int Flags) {
  TString opt = Opt;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("#----------------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
    printf("#   I nsh   SID   Flags  Pln:Pnl:Lay:Str      X        Y        Z       Time     TCorr     eDep   End");
    printf("  DrTime  PrTime  TRes    WDist     WRes        PDG     PDG(M) GenID simID       p        pz\n");
    printf("#----------------------------------------------------------------------------------------------------");
    printf("--------------------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;

  const SimParticle * sim (0);

  int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
  double   mc_mom(-1.);
  double   mc_mom_z(-1.);

  GenId gen_id;

  if (Step) {
    art::Ptr<SimParticle> const& simptr = Step->simParticle();
    art::Ptr<SimParticle> mother        = simptr;

    while(mother->hasParent()) mother = mother->parent();

    sim           = mother.operator ->();

    pdg_id        = simptr->pdgId();
    mother_pdg_id = sim->pdgId();

    if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
    else                         generator_id = -1;

    sim_id        = simptr->id().asInt();
    mc_mom        = Step->momvec().mag();
    mc_mom_z      = Step->momvec().z();
  }

  if ((opt == "") || (opt.Index("data") >= 0)) {
    if (IHit  >= 0) printf("%5i " ,IHit);
    else            printf("      ");

    printf("%3i ",Hit->nStrawHits());

    printf("%5u",Hit->strawId().asUint16());

    if (Flags >= 0) printf(" %08x",Flags);
    else            printf("        ");

    printf(" %3i %3i %3i %3i  %8.3f %8.3f %9.3f %8.3f %8.3f %8.5f  %3i %7.2f %7.2f %5.2f %8.3f %8.3f %10i %10i %5i %5i %8.3f %8.3f\n",
           Hit->strawId().plane(),
           Hit->strawId().panel(),
           Hit->strawId().layer(),
           Hit->strawId().straw(),
           Hit->pos().x(),Hit->pos().y(),Hit->pos().z(),
           Hit->time(),
           Hit->correctedTime(),
           Hit->energyDep(),

           (int) Hit->driftEnd(),
           Hit->driftTime(),
           Hit->propTime(),
           Hit->transRes(),
           Hit->wireDist(),
           Hit->wireRes(),

           pdg_id,
           mother_pdg_id,
           generator_id,
           sim_id,
           mc_mom,
           mc_mom_z);
  }
}
// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(DeltaFinderAna)

}
