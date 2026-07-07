//////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
// this module doesn't do reconstruction
// on input, it takes a list  of flagged hits flags and evaluates the
// performance of the delta electron tagging
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
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
// conditions
#include "Offline/TrackerGeom/inc/Tracker.hh"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

#include "Offline/CalPatRec/inc/HlPrint.hh"
// #include "Offline/CalPatRec/inc/McPart_t.hh"  .. so far, duplicating...

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class DeltaFinderAna : public art::EDAnalyzer {

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>    shCollTag              {Name("shCollTag"             ), Comment("SComboHit collection tag"         ) };
      fhicl::Atom<art::InputTag>    chCollTag              {Name("chCollTag"             ), Comment("ComboHit collection tag"          ) };
      fhicl::Atom<art::InputTag>    sschCollTag            {Name("sschCollTag"           ), Comment("SS ComboHit collection tag"       ) };
      fhicl::Atom<art::InputTag>    sdmcCollTag            {Name("sdmcCollTag"           ), Comment("StrawDigiMC collection Name"      ) };
      fhicl::Atom<int>              debugLevel             {Name("debugLevel"            ), Comment("debug level"                      ) };
      fhicl::Atom<int>              diagLevel              {Name("diagLevel"             ), Comment("diag level"                       ) };
      fhicl::Atom<int>              printElectrons         {Name("printElectrons"        ), Comment("print Electrons"                  ) };
      fhicl::Atom<int>              printElectronsMinNHits {Name("printElectronsMinNHits"), Comment("min N(hits) for printed electrons") };
      fhicl::Atom<float>            printElectronsMaxFReco {Name("printElectronsMaxFReco"), Comment("max F(reco) for printed electrons") };
      fhicl::Atom<int>              printElectronHits      {Name("printElectronHits"     ), Comment("if 1, print electron hits"        ) };
      fhicl::Atom<int>              printComboHits         {Name("printComboHits"        ), Comment("if 1, print combo hits"           ) };
      fhicl::Atom<int>              printSingleComboHits   {Name("printSingleComboHits"  ), Comment("if 1, print single straw ComboH"  ) };
    };

  public:
    enum { kNStations      = 20 };
    enum { kNFaces         =  4 };
    enum { kNPanelsPerFace =  3 };
    enum { kc2_cut         = 50 };
    enum { kPerpRes        = 10 };
    enum { kMaxDxy         = 50 };
    enum { kmax_gap        = 1  };

    enum {
      kNComboHitHistSets = 100,
      kNEventHistSets    =  10,
      kNMcHistSets       = 200,
      kNStrawHitHistSets = 100
    };

    enum { kProtonOrDeut  = 0,
           kLoMomElectron = 1,
           kMdMomElectron = 2,
           kCeMomElectron = 3,
           kHiMomElectron = 4,
           kPositron      = 5,
           kMuon          = 6,
           kOther         = 7,
    };

    struct ComboHitHist_t {
      TH1F*  fNSh;                     // number of straw hits per CH
      TH1F*  fEDep;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNSecondHits;
      TH1F*  fNSecondHitsT;
      TH1F*  fNMc;
      TH1F*  fNHitsDeltaT;
      TH1F*  fNHitsDeltaR;
      TH1F*  fNSSh;                    // total number of straw hits
      TH1F*  fNCh;                     // number of combo hits
      TH1F*  fNSh;                     // number of straw hits in combo hits
      TH1F*  fNChProton;               // N(combo hits) produced by protons
      TH1F*  fNChProtonFP;             // N(proton combo hits flagged as such)
      TH1F*  fNChNonProtonFP;          // N(non-proton combo hits flagged as such)
    };

    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom[2];
      TH1F*  fNCHits;
      TH1F*  fNSSCHits;
      TH2F*  fNCHitsVsMom;
      TH2F*  fNCHitsRVsMom;
      TH2F*  fNCHitsRVsNCHits;
      TH1F*  fNChFlaggedDelta;
      TH1F*  fNChFlaggedProton;
      TH1F*  fFractReco;
      TH2F*  fFractRecoVsNCHits;
    };

    struct StrawHitHist_t {
      TH1F*  fType;
      TH1F*  fTime;
      TH1F*  fMom[2];                       // momentum of the particle which produced the hit
      TH1F*  fEnergyDep;
      TH1F*  fDeltaT;
      TH1F*  fPDGCode;
      TH1F*  fDeltaFlag;
    };

    struct Hist_t {
      ComboHitHist_t* fComboHit[kNComboHitHistSets];
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
      int   fNShFlaggedDelta;     // number of single straw combo hits (SSCH) tagged as delta
      int   fNShFlaggedProton;    // number of single straw combo hits (SSCH) tagged as delta
      int   fTime;                // lowest out of the hit times
      int   fNChFlaggedDelta;
      int   fNChFlaggedProton;
      int   fType;
//-----------------------------------------------------------------------------
// a hit and its flag here could be inconsistent - the flag comes from the 'combo' combohit
//-----------------------------------------------------------------------------
      const SimParticle*               fSim;

      std::vector<const ComboHit*>     fListOfSsComboHits;     // 1-straw combo hits (SSCH)
      //      std::vector<const StrawHitFlag*> fListOfSsComboHitFlags; // a parallel list, flags are those of SS ComboHits

      std::vector<const ComboHit*>     fListOfComboHits;       // "combo" combo hits
      // std::vector<const StrawHitFlag*> fListOfChFlags;         // a parallel list, the flags are those of "Combo" ComboHits

      McPart_t(const SimParticle* Sim = NULL) {
        fSim              = Sim;
        fFirstStation     = 999;
        fLastStation      = -1;
        fNChFlaggedDelta  = 0;
        fNShFlaggedDelta  = 0;
        fNChFlaggedProton = 0;
        fNShFlaggedProton = 0;
        fTime             = 1.e6;
        fType             = MC_Type();
      }

      ~McPart_t() {
      }

      const ComboHit*     SSCh (int I) const { return fListOfSsComboHits[I] ; }
      const ComboHit*     Ch   (int I) const { return fListOfComboHits[I] ; }

      int                 nSsComboHits    () const { return fListOfSsComboHits.size(); }
      int                 nComboHits      () const { return fListOfComboHits.size(); }

      int                 nShFlaggedDelta () const { return fNShFlaggedDelta; }
      int                 nChFlaggedDelta () const { return fNChFlaggedDelta; }

      int                 nShFlaggedProton() const { return fNShFlaggedProton; }
      int                 nChFlaggedProton() const { return fNChFlaggedProton; }

      const StrawHitFlag* ChFlag     (int I) const { return &fListOfComboHits[I]->flag();  }

      float Momentum() const {
        float px = fSim->startMomentum().px();
        float py = fSim->startMomentum().py();
        float pz = fSim->startMomentum().pz();
        return sqrt(px*px+py*py+pz*pz);
      }

      float Time() const { return fTime; }
      float Type() const { return fType; }

                                        // type to index histograms
      int MC_Type() {
        int mc_type(-1);

        float mom = Momentum();

        if      ((fSim->pdgId() == 2212) or (fSim->pdgId() == 1000010020)) mc_type = kProtonOrDeut;
        else if (fSim->pdgId() == 11) {
          if      (mom <  20) mc_type = kLoMomElectron;
          else if (mom <  80) mc_type = kMdMomElectron;
          else if (mom < 110) mc_type = kCeMomElectron;
          else                mc_type = kHiMomElectron;
        }
        else if (fSim->pdgId() == -11)                            mc_type = kPositron;
        else if ((fSim->pdgId() == 13) or (fSim->pdgId() == -13)) mc_type = kMuon;
        else                                                      mc_type = kOther;

        return mc_type;
      }

    };

    struct McHitInfo_t {
      McPart_t*            fMc;
      const  StrawHitFlag* fFlag;
      int                  fType;   // 0:p, 1:ele p<20, 2:ele 20<p<80  3:ele 100<p<110 4:everything else
    };
//-----------------------------------------------------------------------------
// NStations stations, 4-1=3 faces (for hit w/ lower z), 3 panels (for hit w/ lower z)
// 2017-07-27 P.Murat: the 2nd dimension should be 3, right?
// _listOfMcParticles : list of MC particles with at least one digitized tracker hit
// SSH : single-straw hit
//-----------------------------------------------------------------------------
    std::vector<McPart_t*>    _listOfMcParticles;
    std::vector<McHitInfo_t>  _listOfMcHitInfo ; // for each 1-straw hit, pointer to the MC info
    std::vector<McPart_t*>    _listOfProtons;
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag                  _shCollTag;              // straw hits        by "makeSH"
    art::InputTag                  _chCollTag;
    art::InputTag                  _sschCollTag;            // single straw combohits by "makeSH"
    art::InputTag                  _sdmcCollTag;
    int                            _debugLevel;
    int                            _diagLevel;
    int                            _printElectrons;         //
    int                            _printElectronsMinNHits;
    float                          _printElectronsMaxFReco;
    int                            _printElectronHits;      // these are SCHits
    int                            _printComboHits;
    int                            _printSingleComboHits;
//-----------------------------------------------------------------------------
// cache of event or geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl;
    const ComboHitCollection*      _sschColl; // one straw hit per combo hit
    const StrawHitCollection*      _shColl;
    const StrawDigiMCCollection*   _sdmcColl;

    const Tracker*                 _tracker;
    int                            _eventNum;
    int                            _nSingleSH;          // true number of straw hits
    int                            _nComboHits;
                                                        // not the total number of straw hits,
    int                            _nStrawHits;         // but the number of straw hits from combohitse


    int                            _nChFlaggedDelta;
    int                            _nChFlaggedProton;

    int                            _nChDelta;           // total number of CH by e- and e+ P<20 MeV/c
    int                            _nShDelta;
    int                            _nChDeltaFD;         //
    int                            _nChNonDeltaFD;      //

    int                            _nChProton;          // total number of CH by protons
    int                            _nShProton;
    int                            _nChProtonFP;        // number of proton hits flagged as such
    int                            _nChNonProtonFP;     // number of non-proton htis

    HlPrint*                       _hlp;

    StrawHitFlag                   _deltaMask;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit  DeltaFinderAna(const art::EDAnalyzer::Table<Config>& config);
    virtual   ~DeltaFinderAna();

    void      bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookEventHistograms   (EventHist_t*    Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookMcHistograms      (McHist_t*       Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookHistograms();

    void      debug();

    void      fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Hit, McHitInfo_t* McHitInfo);
    void      fillEventHistograms   (EventHist_t*    Hist);
    void      fillMcHistograms      (McHist_t*       Hist, McPart_t* Mc );
    void      fillStrawHitHistograms(StrawHitHist_t* Hist, const StrawHit* Hit, McHitInfo_t* McHitInfo);
    void      fillHistograms();

    bool      findData     (const art::Event&  Evt);
    McPart_t* findParticle (const SimParticle* Sim);
    int       initMcDiag      ();
//-----------------------------------------------------------------------------
// overloaded methods of the base class
//-----------------------------------------------------------------------------
    virtual void beginJob();
    virtual void beginRun(const art::Run&   r);
    virtual void analyze (const art::Event& e);
  };

//-----------------------------------------------------------------------------
  DeltaFinderAna::DeltaFinderAna(const art::EDAnalyzer::Table<Config>& config):
    art::EDAnalyzer(config),
    _shCollTag             (config().shCollTag  ()      ),
    _chCollTag             (config().chCollTag  ()      ),
    _sschCollTag           (config().sschCollTag()      ),
    _sdmcCollTag           (config().sdmcCollTag()      ),
    _debugLevel            (config().debugLevel ()      ),
    _diagLevel             (config().diagLevel  ()      ),
    _printElectrons        (config().printElectrons()   ),
    _printElectronsMinNHits(config().printElectronsMinNHits()),
    _printElectronsMaxFReco(config().printElectronsMaxFReco()),
    _printElectronHits     (config().printElectronHits     ()),
    _printComboHits        (config().printComboHits        ()),
    _printSingleComboHits  (config().printSingleComboHits  ())
  {
    _hlp = HlPrint::Instance();
  }

  DeltaFinderAna::~DeltaFinderAna() {
  }


//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {
    Hist->fNSh      = Dir->make<TH1F>("nsh"  , "N(SH)"        ,   10, 0., 10.);
    Hist->fEDep     = Dir->make<TH1F>("edep" , "edep"         ,  200, 0., 2e-2);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookEventHistograms(EventHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>(Form("event_%02i", HistSet), "Event Number", 100, 0., 100000.);

    Hist->fNSecondHits     = Dir->make<TH1F>(Form("nhit2_%02i" ,HistSet), "N(second hits)", 100, 0., 100.);
    Hist->fNSecondHitsT    = Dir->make<TH1F>(Form("nhit2t_%02i",HistSet), "N(second hits) w/ time selection", 100, 0., 100.);

    Hist->fNMc             = Dir->make<TH1F>(Form("nmc_%02i"       ,HistSet), "N(MC particles)", 100, 0., 1000.);
    Hist->fNHitsDeltaT     = Dir->make<TH1F>(Form("n_delta_ht_%02i",HistSet), "N(delta hits T)", 500, 0., 5000.);
    Hist->fNHitsDeltaR     = Dir->make<TH1F>(Form("n_delta_hr_%02i",HistSet), "N(delta hits R)", 500, 0., 5000.);
    Hist->fNCh             = Dir->make<TH1F>(Form("nch_%02i"       ,HistSet), "N(combo hits)"  ,1000, 0., 10000.);
    Hist->fNSh             = Dir->make<TH1F>(Form("nsh_%02i"       ,HistSet), "N(straw hits)"  ,1000, 0., 10000.);
    Hist->fNSSh            = Dir->make<TH1F>(Form("nssh_%02i"      ,HistSet), "N(1-straw hits)",1000, 0., 10000.);
    Hist->fNChProton       = Dir->make<TH1F>(Form("nch_p_%02i"     ,HistSet), "N(CH prot)"     , 200, 0.,  1000.);
    Hist->fNChProtonFP     = Dir->make<TH1F>(Form("nch_p_fp_%02i"  ,HistSet), "N(CH prot FP)"  , 200, 0.,  1000.);
    Hist->fNChNonProtonFP  = Dir->make<TH1F>(Form("nch_np_fp_%02i" ,HistSet), "N(CH NP FP)"    , 200, 0.,  1000.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookMcHistograms(McHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fPDGCode          = Dir->make<TH1F>("pdg"         , "PDG code"            , 500, -250., 250.);
    Hist->fMom[0]           = Dir->make<TH1F>("mom_0"       , "momentum[0]"         , 500, 0., 500.);
    Hist->fMom[1]           = Dir->make<TH1F>("mom_1"       , "momentum[1]"         , 100, 0.,  20.);
    Hist->fNCHits           = Dir->make<TH1F>("nch"         , "N(combo hits)"       , 200, 0., 200.);
    Hist->fNSSCHits         = Dir->make<TH1F>("nssch"       , "N(SS combo hits)"    , 200, 0., 200.);

    Hist->fNChFlaggedDelta  = Dir->make<TH1F>("nch_fd"      , "N(ch flagged delta)" , 200, 0., 200.);
    Hist->fNChFlaggedProton = Dir->make<TH1F>("nch_fp"      , "N(ch flagged Prot )" , 200, 0., 200.);

    Hist->fNCHitsVsMom      = Dir->make<TH2F>("nch_vs_mom" , "N(hits) vs mom"   , 100,0,100, 200, 0., 200.);
    Hist->fNCHitsRVsMom     = Dir->make<TH2F>("nchr_vs_mom", "N(chits R) vs mom", 100,0,100, 200, 0., 200.);
    Hist->fNCHitsRVsNCHits  = Dir->make<TH2F>("nchr_vs_nh" , "N(chits R) vs cNH", 100,0,100, 100, 0., 100.);
    Hist->fFractReco        = Dir->make<TH1F>("fr"         , "NCHR/NCH"         , 100, 0.,   1.);

    Hist->fFractRecoVsNCHits = Dir->make<TH2F>("fr_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fType      = Dir->make<TH1F>("type"  , "Hit type"        ,   10, 0., 10.);
    Hist->fTime      = Dir->make<TH1F>("time"  , "time"            ,  400, 0., 2000.);
    Hist->fMom[0]    = Dir->make<TH1F>("mom_0" , "Momentum[0]"     ,  100, 0.,  20.);
    Hist->fMom[1]    = Dir->make<TH1F>("mom_1" , "Momentum[1]"     ,  400, 0., 400.);
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
// book combo hit histograms
//-----------------------------------------------------------------------------
    int book_combo_hit_histset[kNComboHitHistSets];
    for (int i=0; i<kNComboHitHistSets; i++) book_combo_hit_histset[i] = 0;

    book_combo_hit_histset[  0] = 1;                // all:
    book_combo_hit_histset[  1] = 1;                // all: prot and deut
    book_combo_hit_histset[  2] = 1;                // all: e-: p<20
    book_combo_hit_histset[  3] = 1;                // all: e- 20<p<80
    book_combo_hit_histset[  4] = 1;                // all: e-: 80<p<110
    book_combo_hit_histset[  5] = 1;                // all: e-  p > 110
    book_combo_hit_histset[  6] = 1;                // all: e+
    book_combo_hit_histset[  7] = 1;                // all: mu- and mu+
    book_combo_hit_histset[  8] = 1;                // all: everything else

    book_combo_hit_histset[ 10] = 1;                // flagged delta:
    book_combo_hit_histset[ 11] = 1;                // flagged delta: prot and deut
    book_combo_hit_histset[ 12] = 1;                // flagged delta: e-: p<20
    book_combo_hit_histset[ 13] = 1;                // flagged delta: e- 20<p<80
    book_combo_hit_histset[ 14] = 1;                // flagged delta: e-: 80<p<110
    book_combo_hit_histset[ 15] = 1;                // flagged delta: e-  p > 110
    book_combo_hit_histset[ 16] = 1;                // flagged delta: e+
    book_combo_hit_histset[ 17] = 1;                // flagged delta: mu- and mu+
    book_combo_hit_histset[ 18] = 1;                // flagged delta: everything else

    book_combo_hit_histset[ 20] = 1;                // not flagged delta:
    book_combo_hit_histset[ 21] = 1;                // not flagged delta: prot and deut
    book_combo_hit_histset[ 22] = 1;                // not flagged delta: e-: p<20
    book_combo_hit_histset[ 23] = 1;                // not flagged delta: e- 20<p<80
    book_combo_hit_histset[ 24] = 1;                // not flagged delta: e-: 80<p<110
    book_combo_hit_histset[ 25] = 1;                // not flagged delta: e-  p > 110
    book_combo_hit_histset[ 26] = 1;                // not flagged delta: e+
    book_combo_hit_histset[ 27] = 1;                // not flagged delta: mu- and mu+
    book_combo_hit_histset[ 28] = 1;                // not flagged delta: everything else

    book_combo_hit_histset[ 30] = 1;                // flagged proton:
    book_combo_hit_histset[ 31] = 1;                // flagged proton: prot and deut
    book_combo_hit_histset[ 32] = 1;                // flagged proton: e-: p<20
    book_combo_hit_histset[ 33] = 1;                // flagged proton: e- 20<p<80
    book_combo_hit_histset[ 34] = 1;                // flagged proton: e-: 80<p<110
    book_combo_hit_histset[ 35] = 1;                // flagged proton: e-  p > 110
    book_combo_hit_histset[ 36] = 1;                // flagged proton: e+
    book_combo_hit_histset[ 37] = 1;                // flagged proton: mu- and mu+
    book_combo_hit_histset[ 38] = 1;                // flagged proton: everything else

    book_combo_hit_histset[ 40] = 1;                // flagged delta or proton:
    book_combo_hit_histset[ 41] = 1;                // flagged delta or proton: prot and deut
    book_combo_hit_histset[ 42] = 1;                // flagged delta or proton: e-: p<20
    book_combo_hit_histset[ 43] = 1;                // flagged delta or proton: e- 20<p<80
    book_combo_hit_histset[ 44] = 1;                // flagged delta or proton: e-: 80<p<110
    book_combo_hit_histset[ 45] = 1;                // flagged delta or proton: e-  p > 110
    book_combo_hit_histset[ 46] = 1;                // flagged delta or proton: e+
    book_combo_hit_histset[ 47] = 1;                // flagged delta or proton: mu- and mu+
    book_combo_hit_histset[ 48] = 1;                // flagged delta or proton: everything else

    for (int i=0; i<kNComboHitHistSets; i++) {
      if (book_combo_hit_histset[i] != 0) {
        sprintf(folder_name,"ch_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fComboHit[i] = new ComboHitHist_t;
        bookComboHitHistograms(_hist.fComboHit[i],i,&tfdir);
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

    book_straw_hit_histset[ 10] = 1;                // flagged delta:
    book_straw_hit_histset[ 11] = 1;                // flagged delta: prot and deut
    book_straw_hit_histset[ 12] = 1;                // flagged delta: e-: p<20
    book_straw_hit_histset[ 13] = 1;                // flagged delta: e- 20<p<80
    book_straw_hit_histset[ 14] = 1;                // flagged delta: e-: 80<p<110
    book_straw_hit_histset[ 15] = 1;                // flagged delta: e-  p > 110
    book_straw_hit_histset[ 16] = 1;                // flagged delta: e+
    book_straw_hit_histset[ 17] = 1;                // flagged delta: mu- and mu+
    book_straw_hit_histset[ 18] = 1;                // flagged delta: everything else

    book_straw_hit_histset[ 20] = 1;                // not flagged delta:
    book_straw_hit_histset[ 21] = 1;                // not flagged delta: prot and deut
    book_straw_hit_histset[ 22] = 1;                // not flagged delta: e-: p<20
    book_straw_hit_histset[ 23] = 1;                // not flagged delta: e- 20<p<80
    book_straw_hit_histset[ 24] = 1;                // not flagged delta: e-: 80<p<110
    book_straw_hit_histset[ 25] = 1;                // not flagged delta: e-  p > 110
    book_straw_hit_histset[ 26] = 1;                // not flagged delta: e+
    book_straw_hit_histset[ 27] = 1;                // not flagged delta: mu- and mu+
    book_straw_hit_histset[ 28] = 1;                // not flagged delta: everything else

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

    book_mc_histset[  1] = 1;                // protons/deuterons
    book_mc_histset[  2] = 1;                // e- p<20
    book_mc_histset[  3] = 1;                // e- 20<p<80
    book_mc_histset[  4] = 1;                // e- 80<p<110
    book_mc_histset[  5] = 1;                // e- p>110      // , 5+ hits, and reco delta
    book_mc_histset[  6] = 1;                // e+       electrons 5+ hits and p < 20
    book_mc_histset[  7] = 1;                // muons
    book_mc_histset[  8] = 1;                // else

    book_mc_histset[ 10] = 1;                // all 5+ hits and reco delta
    book_mc_histset[ 11] = 1;                // protons/deuterons 5+ hits and reco delta
    book_mc_histset[ 12] = 1;                // e- p<20           5+ hits and reco delta
    book_mc_histset[ 13] = 1;                // e- 20<p<80
    book_mc_histset[ 14] = 1;                // e- 80<p<110
    book_mc_histset[ 15] = 1;                // e- p>110      // , 5+ hits, and reco delta
    book_mc_histset[ 16] = 1;                // e+       electrons 5+ hits and p < 20
    book_mc_histset[ 17] = 1;                // muons
    book_mc_histset[ 18] = 1;                // else

    book_mc_histset[ 22] = 1;                // e-    p < 2
    book_mc_histset[ 23] = 1;                // e-    2 < p < 20

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
  void DeltaFinderAna::beginRun(const art::Run& R) {
    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    _tracker = ttHandle.get();
  }

//-----------------------------------------------------------------------------
// assume that, for a ComboHit,all hits are from the same particle
//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Hit, McHitInfo_t* McHitInfo) {
    Hist->fNSh->Fill(Hit->nStrawHits());
    Hist->fEDep->Fill(Hit->energyDep());
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillEventHistograms(EventHist_t* Hist) {
    Hist->fEventNumber->Fill(_eventNum);

//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _listOfMcParticles.size();

    Hist->fNMc->Fill(nmc);
    Hist->fNHitsDeltaT->Fill(_nChDelta);
    Hist->fNHitsDeltaR->Fill(_nChFlaggedDelta);
    Hist->fNChProton->Fill(_nChProton);
    Hist->fNChProtonFP->Fill(_nChProtonFP);
    Hist->fNChNonProtonFP->Fill(_nChNonProtonFP);
    Hist->fNCh->Fill (_nComboHits);
    Hist->fNSh->Fill (_nStrawHits);
    Hist->fNSSh->Fill(_nSingleSH );
  }

//-----------------------------------------------------------------------------
// assume that, for a ComboHit,all hits are from the same particle
//-----------------------------------------------------------------------------
  void  DeltaFinderAna::fillStrawHitHistograms(StrawHitHist_t* Hist, const StrawHit* Hit, McHitInfo_t* McHitInfo) {

    const McPart_t* mc = McHitInfo->fMc;

    Hist->fType->Fill(McHitInfo->fType);
    Hist->fTime->Fill(Hit->time());
    Hist->fMom[0]->Fill(mc->Momentum());
    Hist->fMom[1]->Fill(mc->Momentum());
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
    Hist->fMom[0]->Fill(mom);
    Hist->fMom[1]->Fill(mom);
    Hist->fNCHits->Fill(Mc->nComboHits());
    Hist->fNSSCHits->Fill(Mc->nSsComboHits());
    Hist->fNCHitsVsMom->Fill(mom,Mc->nComboHits());

    Hist->fNChFlaggedDelta->Fill (Mc->nChFlaggedDelta());
    Hist->fNChFlaggedProton->Fill(Mc->nChFlaggedProton());

    Hist->fNCHitsRVsMom->Fill(mom,Mc->nChFlaggedDelta());
    Hist->fNCHitsRVsNCHits->Fill(Mc->nComboHits(),Mc->nChFlaggedDelta());

    float freco = Mc->nChFlaggedDelta()/(Mc->nComboHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNCHits->Fill(Mc->nComboHits(),freco);
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
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch  = &_chColl ->at(i);
      //      const StrawHitFlag* chf = &_chfColl->at(i);
      const StrawHitFlag* chf = &ch->flag();

      //      int nsh = ch->nStrawHits();
//-----------------------------------------------------------------------------
// 'mark' CH by its first straw hit - with some unavoidable uncertainties
//-----------------------------------------------------------------------------
      int ind                  = ch->index(0);
      // const StrawHit*     sh   = &_shColl->at(ind);
      const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = sgs->simParticle().get();

      McPart_t mc(sim);
      McHitInfo_t    mc_hit_info;

      int mc_type = mc.Type();

      mc_hit_info.fMc   = &mc;
      mc_hit_info.fType = mc_type;
      //      mc_hit_info.fFlag = &_chfColl->at(i);
      mc_hit_info.fFlag = chf;

      fillComboHitHistograms(_hist.fComboHit[0],ch,&mc_hit_info);  // all

      if      (mc_type == kProtonOrDeut ) fillComboHitHistograms(_hist.fComboHit[1],ch,&mc_hit_info);
      else if (mc_type == kLoMomElectron) fillComboHitHistograms(_hist.fComboHit[2],ch,&mc_hit_info);
      else if (mc_type == kMdMomElectron) fillComboHitHistograms(_hist.fComboHit[3],ch,&mc_hit_info);
      else if (mc_type == kCeMomElectron) fillComboHitHistograms(_hist.fComboHit[4],ch,&mc_hit_info);
      else if (mc_type == kHiMomElectron) fillComboHitHistograms(_hist.fComboHit[5],ch,&mc_hit_info);
      else if (mc_type == kPositron     ) fillComboHitHistograms(_hist.fComboHit[6],ch,&mc_hit_info);
      else if (mc_type == kMuon         ) fillComboHitHistograms(_hist.fComboHit[7],ch,&mc_hit_info);
      else                                fillComboHitHistograms(_hist.fComboHit[8],ch,&mc_hit_info);

      int delta  = chf->hasAnyProperty(StrawHitFlag::bkg);
      int proton = (chf->hasAnyProperty(StrawHitFlag::energysel) == false);

      if (delta) {
        fillComboHitHistograms(_hist.fComboHit[10],ch,&mc_hit_info);  // all

        if      (mc_type == kProtonOrDeut ) fillComboHitHistograms(_hist.fComboHit[11],ch,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillComboHitHistograms(_hist.fComboHit[12],ch,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillComboHitHistograms(_hist.fComboHit[13],ch,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillComboHitHistograms(_hist.fComboHit[14],ch,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillComboHitHistograms(_hist.fComboHit[15],ch,&mc_hit_info);
        else if (mc_type == kPositron     ) fillComboHitHistograms(_hist.fComboHit[16],ch,&mc_hit_info);
        else if (mc_type == kMuon         ) fillComboHitHistograms(_hist.fComboHit[17],ch,&mc_hit_info);
        else                                fillComboHitHistograms(_hist.fComboHit[18],ch,&mc_hit_info);
      }
      else {
        fillComboHitHistograms(_hist.fComboHit[20],ch,&mc_hit_info);  // all

        if      (mc_type == kProtonOrDeut ) fillComboHitHistograms(_hist.fComboHit[21],ch,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillComboHitHistograms(_hist.fComboHit[22],ch,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillComboHitHistograms(_hist.fComboHit[23],ch,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillComboHitHistograms(_hist.fComboHit[24],ch,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillComboHitHistograms(_hist.fComboHit[25],ch,&mc_hit_info);
        else if (mc_type == kPositron     ) fillComboHitHistograms(_hist.fComboHit[26],ch,&mc_hit_info);
        else if (mc_type == kMuon         ) fillComboHitHistograms(_hist.fComboHit[27],ch,&mc_hit_info);
        else                                fillComboHitHistograms(_hist.fComboHit[28],ch,&mc_hit_info);
      }

      if (proton) {
        fillComboHitHistograms(_hist.fComboHit[30],ch,&mc_hit_info);  // all

        if      (mc_type == kProtonOrDeut ) fillComboHitHistograms(_hist.fComboHit[31],ch,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillComboHitHistograms(_hist.fComboHit[32],ch,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillComboHitHistograms(_hist.fComboHit[33],ch,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillComboHitHistograms(_hist.fComboHit[34],ch,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillComboHitHistograms(_hist.fComboHit[35],ch,&mc_hit_info);
        else if (mc_type == kPositron     ) fillComboHitHistograms(_hist.fComboHit[36],ch,&mc_hit_info);
        else if (mc_type == kMuon         ) fillComboHitHistograms(_hist.fComboHit[37],ch,&mc_hit_info);
        else                                fillComboHitHistograms(_hist.fComboHit[38],ch,&mc_hit_info);
      }

      if (delta or proton) {
        fillComboHitHistograms(_hist.fComboHit[40],ch,&mc_hit_info);  // all

        if      (mc_type == kProtonOrDeut ) fillComboHitHistograms(_hist.fComboHit[41],ch,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillComboHitHistograms(_hist.fComboHit[42],ch,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillComboHitHistograms(_hist.fComboHit[43],ch,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillComboHitHistograms(_hist.fComboHit[44],ch,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillComboHitHistograms(_hist.fComboHit[45],ch,&mc_hit_info);
        else if (mc_type == kPositron     ) fillComboHitHistograms(_hist.fComboHit[46],ch,&mc_hit_info);
        else if (mc_type == kMuon         ) fillComboHitHistograms(_hist.fComboHit[47],ch,&mc_hit_info);
        else                                fillComboHitHistograms(_hist.fComboHit[48],ch,&mc_hit_info);
      }
   }

    for (int i=0; i<_nSingleSH; i++) {
      const ComboHit*     ssch  = &_sschColl ->at(i);
      //      const StrawHitFlag* sschf = &_sschfColl->at(i);

      int ind = ssch->index(0);

      const StrawHit*     sh   = &_shColl->at(ind);
      const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = sgs->simParticle().get();

      McPart_t mc(sim);
      McHitInfo_t    mc_hit_info;

      int mc_type = mc.Type();

      mc_hit_info.fMc   = &mc;
      mc_hit_info.fType = mc_type;
      mc_hit_info.fFlag = &ssch->flag();
//-----------------------------------------------------------------------------
// fill single straw histograms
//-----------------------------------------------------------------------------
      fillStrawHitHistograms(_hist.fStrawHit[0],sh,&mc_hit_info);  // all

      if      (mc_type == kProtonOrDeut ) fillStrawHitHistograms(_hist.fStrawHit[1],sh,&mc_hit_info);
      else if (mc_type == kLoMomElectron) fillStrawHitHistograms(_hist.fStrawHit[2],sh,&mc_hit_info);
      else if (mc_type == kMdMomElectron) fillStrawHitHistograms(_hist.fStrawHit[3],sh,&mc_hit_info);
      else if (mc_type == kCeMomElectron) fillStrawHitHistograms(_hist.fStrawHit[4],sh,&mc_hit_info);
      else if (mc_type == kHiMomElectron) fillStrawHitHistograms(_hist.fStrawHit[5],sh,&mc_hit_info);
      else if (mc_type == kPositron     ) fillStrawHitHistograms(_hist.fStrawHit[6],sh,&mc_hit_info);
      else if (mc_type == kMuon         ) fillStrawHitHistograms(_hist.fStrawHit[7],sh,&mc_hit_info);
      else                                fillStrawHitHistograms(_hist.fStrawHit[8],sh,&mc_hit_info);

      const StrawHitFlag* flag = mc_hit_info.fFlag;

      int delta  = flag->hasAllProperties(StrawHitFlag::bkg);

      if (delta) {
//-----------------------------------------------------------------------------
// set 2: all hits flagged as delta electrons
//-----------------------------------------------------------------------------
        fillStrawHitHistograms(_hist.fStrawHit[10],sh,&mc_hit_info);

        if      (mc_type == kProtonOrDeut ) fillStrawHitHistograms(_hist.fStrawHit[11],sh,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillStrawHitHistograms(_hist.fStrawHit[12],sh,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillStrawHitHistograms(_hist.fStrawHit[13],sh,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillStrawHitHistograms(_hist.fStrawHit[14],sh,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillStrawHitHistograms(_hist.fStrawHit[15],sh,&mc_hit_info);
        else if (mc_type == kPositron     ) fillStrawHitHistograms(_hist.fStrawHit[16],sh,&mc_hit_info);
        else if (mc_type == kMuon         ) fillStrawHitHistograms(_hist.fStrawHit[17],sh,&mc_hit_info);
        else                                fillStrawHitHistograms(_hist.fStrawHit[18],sh,&mc_hit_info);
      }
      else {
        fillStrawHitHistograms(_hist.fStrawHit[20],sh,&mc_hit_info);

        if      (mc_type == kProtonOrDeut ) fillStrawHitHistograms(_hist.fStrawHit[21],sh,&mc_hit_info);
        else if (mc_type == kLoMomElectron) fillStrawHitHistograms(_hist.fStrawHit[22],sh,&mc_hit_info);
        else if (mc_type == kMdMomElectron) fillStrawHitHistograms(_hist.fStrawHit[23],sh,&mc_hit_info);
        else if (mc_type == kCeMomElectron) fillStrawHitHistograms(_hist.fStrawHit[24],sh,&mc_hit_info);
        else if (mc_type == kHiMomElectron) fillStrawHitHistograms(_hist.fStrawHit[25],sh,&mc_hit_info);
        else if (mc_type == kPositron     ) fillStrawHitHistograms(_hist.fStrawHit[26],sh,&mc_hit_info);
        else if (mc_type == kMuon         ) fillStrawHitHistograms(_hist.fStrawHit[27],sh,&mc_hit_info);
        else                                fillStrawHitHistograms(_hist.fStrawHit[28],sh,&mc_hit_info);
      }

      //      loc++;
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _listOfMcParticles.size();

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _listOfMcParticles.at(i);
      const SimParticle* sim = mc->fSim;

      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if      (mc->fType == kProtonOrDeut ) fillMcHistograms(_hist.fMc[1],mc);
      else if (mc->fType == kLoMomElectron) fillMcHistograms(_hist.fMc[2],mc);
      else if (mc->fType == kMdMomElectron) fillMcHistograms(_hist.fMc[3],mc);
      else if (mc->fType == kCeMomElectron) fillMcHistograms(_hist.fMc[4],mc);
      else if (mc->fType == kHiMomElectron) fillMcHistograms(_hist.fMc[5],mc);
      else if (mc->fType == kPositron     ) fillMcHistograms(_hist.fMc[6],mc);
      else if (mc->fType == kMuon         ) fillMcHistograms(_hist.fMc[7],mc);
      else                                  fillMcHistograms(_hist.fMc[8],mc);

      if (mc->nComboHits() >= 5) {
        if      (mc->fType == kProtonOrDeut ) fillMcHistograms(_hist.fMc[11],mc);
        else if (mc->fType == kLoMomElectron) fillMcHistograms(_hist.fMc[12],mc);
        else if (mc->fType == kMdMomElectron) fillMcHistograms(_hist.fMc[13],mc);
        else if (mc->fType == kCeMomElectron) fillMcHistograms(_hist.fMc[14],mc);
        else if (mc->fType == kHiMomElectron) fillMcHistograms(_hist.fMc[15],mc);
        else if (mc->fType == kPositron     ) fillMcHistograms(_hist.fMc[16],mc);
        else if (mc->fType == kMuon         ) fillMcHistograms(_hist.fMc[17],mc);
        else                                  fillMcHistograms(_hist.fMc[18],mc);
      }

      if (mc->fType == kLoMomElectron) {
        float mom = sim->startMomXYZT().P();
        if (mom < 1.5) fillMcHistograms(_hist.fMc[22],mc);
        else           fillMcHistograms(_hist.fMc[23],mc);
      }
//-----------------------------------------------------------------------------
// a closer look at misreconstructed delta electrons
//-----------------------------------------------------------------------------
      float fr_ch = mc->nChFlaggedDelta()/(mc->nComboHits()+1.e-3);

      if ((mc->Momentum() < 5) && (mc->Time() > 550) && (mc->nComboHits() > 40) && (fr_ch < 0.5)) {
        printf(" event: %6i missed delta: sim.id = %10li mom = %10.3f time= %9.3f ",
               _eventNum, sim->id().asInt(), mc->Momentum(), mc->Time());
        printf("nch:nch(delta) = %3i:%3i stations(first:last): %2i:%2i",
               mc->nComboHits(), mc->nChFlaggedDelta(), mc->fFirstStation, mc->fLastStation);
        printf(" fraction(ch): %6.3f\n",fr_ch);
      }
    }
  }

//-----------------------------------------------------------------------------
  DeltaFinderAna::McPart_t* DeltaFinderAna::findParticle(const SimParticle* Sim) {
    McPart_t* found(0);

    int n = _listOfMcParticles.size();

    for (int i=0; i<n; i++) {
      McPart_t* mc = _listOfMcParticles.at(i);
      if (mc->fSim == Sim) {
        found = mc;
        break;
      }
    }

    return found;
  }

//-----------------------------------------------------------------------------
bool DeltaFinderAna::findData(const art::Event& Evt) {
    _chColl     = nullptr;
    _sdmcColl   = nullptr;
    // _chfColl    = nullptr;

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();
    _nComboHits = _chColl->size();

    auto shcH   = Evt.getValidHandle<StrawHitCollection>(_shCollTag);
    _shColl     = shcH.product();
    _nSingleSH  = _shColl->size();

    auto sschcH = Evt.getValidHandle<ComboHitCollection>(_sschCollTag);
    _sschColl   = sschcH.product();

    // auto chfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_chfCollTag);
    // _chfColl    = chfcH.product();

    // auto shfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_shfCollTag);
    // _shfColl    = shfcH.product();

    auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
    _sdmcColl   = sdmccH.product();

    //    return (_chColl != 0) && (_chfColl != 0) && (_shfColl != 0) && (_sdmcColl != 0) ;
    return (_chColl != 0) && (_sdmcColl != 0) ;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker
// it should be shorter than the list of all MC particles
//-----------------------------------------------------------------------------
  int DeltaFinderAna::initMcDiag() {

    int n = _listOfMcParticles.size();
    for (int i=0; i<n; i++) {
      McPart_t* p = _listOfMcParticles.at(i);
      delete p;
    }

    _listOfMcParticles.clear();
    _listOfMcHitInfo.clear();
    _listOfProtons.clear();
//-----------------------------------------------------------------------------
// count the number of participating straw hits (remember, it could be less than
// the total number of straw hits :participating are only the straw hits from
// combo hits
//-----------------------------------------------------------------------------
    _nStrawHits = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit* ch = &_chColl->at(i);
      _nStrawHits       += ch->nStrawHits();
    }
//-----------------------------------------------------------------------------
// create list of MC particles - all particles with at least one digitized
// single straw tracker hit
//-----------------------------------------------------------------------------
    for (int i=0; i<_nSingleSH; i++) {
      const StrawDigiMC*  sdmc = &_sdmcColl->at(i);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = sgs->simParticle().get();

      McPart_t* mc = findParticle(sim);

      if (mc == NULL) {
                                        // add new particle
        mc = new McPart_t(sim);
        _listOfMcParticles.push_back(mc);
        if (mc->fType == kProtonOrDeut) _listOfProtons.push_back(mc);
      }
    }
//-----------------------------------------------------------------------------
// count total numbers of delta-electron and proton hits
//-----------------------------------------------------------------------------
    _nShDelta  = 0;
    _nChDelta  = 0;
    _nShProton = 0;
    _nChProton = 0;

    int nmc    = _listOfMcParticles.size();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc    = _listOfMcParticles.at(i);
      if (mc->fType == kLoMomElectron) {
//-----------------------------------------------------------------------------
// call this "a delta electron"
//-----------------------------------------------------------------------------
        _nShDelta += mc->nSsComboHits();
        _nChDelta += mc->nComboHits();
      }
      else if (mc->fType == kProtonOrDeut) {
        _nShProton += mc->nSsComboHits();
        _nChProton += mc->nComboHits();
      }
    }
//-----------------------------------------------------------------------------
// initialize single-straw level information
//-----------------------------------------------------------------------------
    _listOfMcHitInfo.resize(_nSingleSH);
    for (int i=0; i<_nSingleSH; i++) {
      const StrawDigiMC*  sdmc = &_sdmcColl->at(i);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = sgs->simParticle().get();
      McPart_t*           mc   = findParticle(sim);
//-----------------------------------------------------------------------------
// at this point, all MC particles should already be registered
//-----------------------------------------------------------------------------
      const ComboHit* ssch    = &_sschColl->at(i);
      // const StrawHitFlag* shf = &_shfColl->at(i);
      // const StrawHitFlag* shf = &ssch->flag();
//-----------------------------------------------------------------------------
// the SSCH flags from the flag collection could be redefined and different
// from the ones stored in the hit payload
//-----------------------------------------------------------------------------
      mc->fListOfSsComboHits.push_back(ssch);     // these come with their own flags in payload
      //      mc->fListOfSsComboHitFlags.push_back(shf);

      int station = ssch->strawId().getStation();

      if (station < mc->fFirstStation) mc->fFirstStation = station;
      if (station > mc->fLastStation ) mc->fLastStation  = station;

      if (ssch->correctedTime() < mc->fTime) mc->fTime = ssch->correctedTime();

      McHitInfo_t* mc_hit_info = &_listOfMcHitInfo.at(i);

      mc_hit_info->fMc   = mc;
      mc_hit_info->fType = mc->Type();
    }
//-----------------------------------------------------------------------------
// 1) 'shf' is also a mu2e::ComboHit flag, despite the name
// so some straw hits may have flags not consistent with their origin
//-----------------------------------------------------------------------------
    _nChFlaggedDelta  = 0;
    _nChFlaggedProton = 0;

    // int loc = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch   = &_chColl->at(i);
      //      const StrawHitFlag* chf  = &_chfColl->at(i);
      const StrawHitFlag* chf  = &ch->flag();

      int       ind = ch->index(0);
      McPart_t* mc  = _listOfMcHitInfo.at(ind).fMc;

      mc->fListOfComboHits.push_back(ch);
      // mc->fListOfChFlags.push_back  (chf);
//-----------------------------------------------------------------------------
// delta and proton flagging is done using combo hits
//-----------------------------------------------------------------------------
      int flagged_delta  = chf->hasAnyProperty(StrawHitFlag::bkg);
      int flagged_proton = (chf->hasAnyProperty(StrawHitFlag::energysel) == 0);

      if (flagged_delta ) {
        _nChFlaggedDelta      += 1;
        mc->fNChFlaggedDelta  += 1;
      }

      if (flagged_proton) {
        _nChFlaggedProton     += 1;
        mc->fNChFlaggedProton += 1;
      }

      int nsh = ch->nStrawHits();
      for (int ish=0; ish<nsh; ish++) {
        if (flagged_delta ) mc->fNShFlaggedDelta  += 1;
        if (flagged_proton) mc->fNShFlaggedProton += 1;
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::analyze(const art::Event& E) {

    _eventNum = E.event();
    if (_debugLevel) {
      printf("* >>> DeltaFinderAna::%s event number: %10i\n",__func__,_eventNum);
    }

    _hlp->SetEvent(&E);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(E)) {
      throw cet::exception("RECO")
        << "mu2e::DeltaFinderAna_module::produce: missing data" << endl;
    }

    initMcDiag      ();
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
      printf("* DeltaFinderAna::debug : printElectrons\n");
      int nmc = _listOfMcParticles.size();

      for (int i=0; i<nmc; i++) {
        McPart_t* mc = _listOfMcParticles.at(i);
        const SimParticle* sim = mc->fSim;

        if (((sim->pdgId() == PDGCode::e_minus) or (sim->pdgId() == PDGCode::e_plus)) and (mc->nComboHits() >= _printElectronsMinNHits)) {

          float fr_sh = mc->nShFlaggedDelta()/(mc->nSsComboHits()+1.e-3);
          float fr_ch = mc->nChFlaggedDelta()/(mc->nComboHits  ()+1.e-3);

          if (fr_ch < _printElectronsMaxFReco) {

            printf("* electron: sim.id: %6li mom: %8.3f time: %8.2f ",
                   sim->id().asInt(), mc->Momentum(), mc->Time());
            printf("nch:nchd:fr: %3i:%3i:%6.3f  nsh:nshd:fr %3i:%3i:%6.3f stations: %2i:%2i",
                   mc->nComboHits(),mc->nChFlaggedDelta(),fr_ch,
                   mc->nSsComboHits  (),mc->nShFlaggedDelta(),fr_sh,
                   mc->fFirstStation, mc->fLastStation);
            printf(" \n");
//-----------------------------------------------------------------------------
// check if want to print hits - print single-straw combo hits....
//-----------------------------------------------------------------------------
            if (_printElectronHits != 0) {
                                        // these are 'combo' ComboHits
              int nh  = mc->nComboHits();

              _hlp->printComboHit(0,0,"banner",0,0);

              for (int i=0; i<nh; i++) {
                const ComboHit* ch = mc->Ch(i);
                                        // flags here are the combo hit flags, not single-straw hit flags
                int flags  = *((int*) mc->ChFlag(i));
                                        // use sim_id of the first SSCH ...
                int ind    = ch->index(0);

                const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
                const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
                _hlp->printComboHit(ch,sgs,"data",ind,flags);
              }
            }
          }
        }
      }
    }

    if (_printComboHits) {
      printf("* DeltaFinderAna::debug : print ComboHits \n");
//-----------------------------------------------------------------------------
// print ComboHits with flags stored in the external hit collection
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_chCollTag.encode().data(),
                                    "",                            // _chfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }

    if (_printSingleComboHits) {
      printf("* DeltaFinderAna::debug : single straw ComboHits tag:  %s\n",_shCollTag.encode().data());
//-----------------------------------------------------------------------------
// print single-straw ComboHits with flags stored in the external hit collection
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_shCollTag.encode().data(),
                                    "",                            // _shfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }
  }

// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(DeltaFinderAna)

}
