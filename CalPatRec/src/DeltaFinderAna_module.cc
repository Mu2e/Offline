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
#include "art/Framework/Core/EDAnalyzer.h"
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

  class DeltaFinderAna : public art::EDAnalyzer {

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>    shCollTag              {Name("shCollTag"             ), Comment("SComboHit collection tag"         ) };
      fhicl::Atom<art::InputTag>    chCollTag              {Name("chCollTag"             ), Comment("ComboHit collection tag"          ) };
      fhicl::Atom<art::InputTag>    sschCollTag            {Name("sschCollTag"           ), Comment("SS ComboHit collection tag"       ) };
      fhicl::Atom<art::InputTag>    chfCollTag             {Name("chfCollTag"            ), Comment("ComboHit flag collection tag"     ) };
      fhicl::Atom<art::InputTag>    shfCollTag             {Name("shfCollTag"            ), Comment("SSC Hit  flag collection tag"     ) };
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
      kNEventHistSets    =  10,
      kNStrawHitHistSets = 100,
      kNMcHistSets       = 200
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

  protected:

    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom[2];
      TH1F*  fNCHits;
      TH1F*  fNSSCHits;
      TH2F*  fNCHitsVsMom;
      TH2F*  fNCHitsRVsMom;
      TH2F*  fNCHitsRVsNCHits;
      TH1F*  fNCHitsDelta;
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
      int   fNShTaggedDelta;      // number of single straw combo hits (SSCH) tagged as delta
      int   fTime;                // lowest out of the hit times
      int   fNChTaggedDelta;
      int   fType;

      const SimParticle*               fSim;
//-----------------------------------------------------------------------------
// a hit and its flag here could be inconsistent - the flag comes from the 'combo' combohit
//-----------------------------------------------------------------------------
      std::vector<const ComboHit*>     fListOfSSCHits; // 1-straw combo hits (SSCH)
      std::vector<const StrawHitFlag*> fListOfFlags;   // a parallel list, but the flags are those of "Combo" ComboHits

      std::vector<const ComboHit*>     fListOfComboHits; // combo hits (SSCH)
      std::vector<const StrawHitFlag*> fListOfChFlags;   // a parallel list, but the flags are those of "Combo" ComboHits

      McPart_t(const SimParticle* Sim = NULL) {
        fSim          = Sim;
        fFirstStation = 999;
        fLastStation  = -1;
        fNChTaggedDelta = 0;
        fNShTaggedDelta = 0;
        fTime           = 1.e6;
        fType           = MC_Type();
      }

      ~McPart_t() {
      }

      const ComboHit*     SSCh (int I) const { return fListOfSSCHits[I] ;    }
      const ComboHit*     Ch   (int I) const { return fListOfComboHits[I] ;    }

                                        // these are single-hit SCH's
      int                 NSSCHits  () const { return fListOfSSCHits.size(); }
      int                 NComboHits() const { return fListOfComboHits.size(); }

      int                 NShTaggedDelta() const { return fNShTaggedDelta; }
      int                 NChTaggedDelta() const { return fNChTaggedDelta; }

      const StrawHitFlag* ShFlag(int I) const { return fListOfFlags[I];    }
      const StrawHitFlag* ChFlag(int I) const { return fListOfChFlags[I];  }

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
    art::InputTag                  _shCollTag;              // straw hits        by "makeSH"
    art::InputTag                  _chCollTag;
    art::InputTag                  _sschCollTag;            // single straw combohits by "makeSH"
    art::InputTag                  _shfCollTag;             // 1-straw  ComboHit flags
    art::InputTag                  _chfCollTag;             // combined ComboHit flags
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
    const StrawHitFlagCollection*  _chfColl;
    const StrawDigiMCCollection*   _sdmcColl;

    const Tracker*                 _tracker;
    int                            _eventNum;
    int                            _nSingleSH;  // true number of straw hits
    int                            _nComboHits;
    int                            _nStrawHits; // not the total number of straw hits, but the
                                                // number of straw hits from combohits

    int                            fNHitsDeltaTot;
    int                            fNHitsDeltaReco;

    HlPrint*                       _hlp;

    StrawHitFlag                   _deltaMask;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit  DeltaFinderAna(const art::EDAnalyzer::Table<Config>& config);
    virtual   ~DeltaFinderAna();

    int       associateMcTruth();

    void      bookEventHistograms   (EventHist_t*    Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookStrawHitHistograms(StrawHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookMcHistograms      (McHist_t*       Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookHistograms();

    void      debug();

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
    _shfCollTag            (config().shfCollTag ()      ),
    _chfCollTag            (config().chfCollTag ()      ),
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
// form a list of MC particles with hits in the tracker
//-----------------------------------------------------------------------------
  int DeltaFinderAna::associateMcTruth() {
//-----------------------------------------------------------------------------
// for each MC electron calculate the number of reconstructed hits
//-----------------------------------------------------------------------------
    int nmc    = _list_of_mc_particles.size();

    fNHitsDeltaTot  = 0;

    for (int i=0; i<nmc; i++) {
      McPart_t* mc    = _list_of_mc_particles.at(i);

      if (mc->fType == kLoMomElectron) {
//-----------------------------------------------------------------------------
// call this "a delta electron"
//-----------------------------------------------------------------------------
        fNHitsDeltaTot += mc->NSSCHits();
      }
    }

    return 0;
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
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAna::bookMcHistograms(McHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fPDGCode         = Dir->make<TH1F>("pdg"         , "PDG code"        , 500, -250., 250.);
    Hist->fMom[0]          = Dir->make<TH1F>("mom_0"       , "momentum[0]"     , 500, 0., 500.);
    Hist->fMom[1]          = Dir->make<TH1F>("mom_1"       , "momentum[1]"     , 100, 0.,  20.);
    Hist->fNCHits          = Dir->make<TH1F>("nch"         , "N(combo hits)"   , 200, 0., 200.);
    Hist->fNSSCHits        = Dir->make<TH1F>("nssch"       , "N(SS combo hits)", 200, 0., 200.);
    Hist->fNCHitsDelta     = Dir->make<TH1F>("nchitsr"     , "N(chits reco)"   , 200, 0., 200.);

    Hist->fNCHitsVsMom     = Dir->make<TH2F>("nch_vs_mom" , "N(hits) vs mom"   , 100,0,100, 200, 0., 200.);
    Hist->fNCHitsRVsMom    = Dir->make<TH2F>("nchr_vs_mom", "N(chits R) vs mom", 100,0,100, 200, 0., 200.);
    Hist->fNCHitsRVsNCHits = Dir->make<TH2F>("nchr_vs_nh" , "N(chits R) vs cNH", 100,0,100, 100, 0., 100.);
    Hist->fFractReco       = Dir->make<TH1F>("fr"         , "NCHR/NCH"         , 100, 0.,   1.);

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
  void  DeltaFinderAna::fillEventHistograms(EventHist_t* Hist) {
    Hist->fEventNumber->Fill(_eventNum);

//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();

    Hist->fNMc->Fill(nmc);
    Hist->fNHitsDeltaT->Fill(fNHitsDeltaTot);
    Hist->fNHitsDeltaR->Fill(fNHitsDeltaReco);
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
    Hist->fNCHits->Fill(Mc->NComboHits());
    Hist->fNSSCHits->Fill(Mc->NSSCHits());
    Hist->fNCHitsVsMom->Fill(mom,Mc->NComboHits());
    Hist->fNCHitsDelta->Fill(Mc->NChTaggedDelta());
    Hist->fNCHitsRVsMom->Fill(mom,Mc->NChTaggedDelta());
    Hist->fNCHitsRVsNCHits->Fill(Mc->NComboHits(),Mc->NChTaggedDelta());

    float freco = Mc->NChTaggedDelta()/(Mc->NComboHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNCHits->Fill(Mc->NComboHits(),freco);
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

        McPart_t mc(sim);
        McHitInfo_t    mc_hit_info;

        int mc_type = mc.Type();

        mc_hit_info.fMc   = &mc;
        mc_hit_info.fType = mc_type;
        mc_hit_info.fFlag = &_chfColl->at(i);

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

        // int edepOK = flag->hasAllProperties(StrawHitFlag::energysel);
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
      if      (mc->fType == kProtonOrDeut ) fillMcHistograms(_hist.fMc[1],mc);
      else if (mc->fType == kLoMomElectron) fillMcHistograms(_hist.fMc[2],mc);
      else if (mc->fType == kMdMomElectron) fillMcHistograms(_hist.fMc[3],mc);
      else if (mc->fType == kCeMomElectron) fillMcHistograms(_hist.fMc[4],mc);
      else if (mc->fType == kHiMomElectron) fillMcHistograms(_hist.fMc[5],mc);
      else if (mc->fType == kPositron     ) fillMcHistograms(_hist.fMc[6],mc);
      else if (mc->fType == kMuon         ) fillMcHistograms(_hist.fMc[7],mc);
      else                                  fillMcHistograms(_hist.fMc[8],mc);

      if (mc->NComboHits() >= 5) {
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
      float fr_ch = mc->NChTaggedDelta()/(mc->NComboHits()+1.e-3);

      if ((mc->Momentum() < 5) && (mc->Time() > 550) && (mc->NComboHits() > 40) && (fr_ch < 0.5)) {
        printf(" event: %6i missed delta: sim.id = %10li mom = %10.3f time= %9.3f ",
               _eventNum, sim->id().asInt(), mc->Momentum(), mc->Time());
        printf("nch:nch(delta) = %3i:%3i stations(first:last): %2i:%2i",
               mc->NComboHits(), mc->NChTaggedDelta(), mc->fFirstStation, mc->fLastStation);
        printf(" fraction(ch): %6.3f\n",fr_ch);
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
bool DeltaFinderAna::findData(const art::Event& Evt) {
    _chColl     = nullptr;
    _sdmcColl   = nullptr;
    _chfColl    = nullptr;

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();
    _nComboHits = _chColl->size();

    auto shcH   = Evt.getValidHandle<StrawHitCollection>(_shCollTag);
    _shColl     = shcH.product();
    _nSingleSH  = _shColl->size();

    auto sschcH = Evt.getValidHandle<ComboHitCollection>(_sschCollTag);
    _sschColl   = sschcH.product();

    auto chfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_chfCollTag);
    _chfColl    = chfcH.product();

    auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
    _sdmcColl   = sdmccH.product();

    return (_chColl != 0) && (_chfColl != 0) && (_sdmcColl != 0) ;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker
// it should be shorter than the list of all MC particles
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
      const ComboHit* ch = &_chColl->at(i);
      _nStrawHits       += ch->nStrawHits();
    }

    _list_of_mc_hit_info.resize(_nSingleSH);
//-----------------------------------------------------------------------------
// create list of MC particles - of all particles which made at least one hit
//-----------------------------------------------------------------------------
    for (int i=0; i<_nSingleSH; i++) {
      const StrawDigiMC*  sdmc = &_sdmcColl->at(i);
      const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
      const SimParticle*  sim  = sgs->simParticle().get();

      McPart_t* mc = findParticle(sim);

      if (mc == NULL) {
                                        // add new particle
        mc = new McPart_t(sim);
        _list_of_mc_particles.push_back(mc);
      }
    }

    fNHitsDeltaReco = 0;
//-----------------------------------------------------------------------------
// 'shf' is also a mu2e::ComboHit flag, despite the name
// so some straw hits may have flags not consistent with their origin
//-----------------------------------------------------------------------------
    int loc = 0;
    for (int i=0; i<_nComboHits; i++) {
      const ComboHit*     ch   = &_chColl->at(i);
      const StrawHitFlag* chf  = &_chfColl->at(i);

      int nsh = ch->nStrawHits();
      int ch_counted = 0;
      for (int ish=0; ish<nsh; ish++) {
        int ind = ch->indexArray().at(ish);
        const ComboHit* ssch = &_sschColl->at(ind);

        const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
        const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();
        const SimParticle*  sim  = sgs->simParticle().get();
//-----------------------------------------------------------------------------
// at this poit, all MC particles should already be registered
//-----------------------------------------------------------------------------
        McPart_t* mc = findParticle(sim);

        mc->fListOfSSCHits.push_back(ssch);

        int station = ssch->strawId().getStation();

        if (station < mc->fFirstStation) mc->fFirstStation = station;
        if (station > mc->fLastStation ) mc->fLastStation  = station;

        if (ssch->correctedTime() < mc->fTime) mc->fTime = ssch->correctedTime();
//-----------------------------------------------------------------------------
// remember, this is a combohit flag !
// delta-wise, all straw hits corresponding to the same combo hit,
// should have the same flags
//-----------------------------------------------------------------------------
        mc->fListOfFlags.push_back(chf);

        McHitInfo_t* mc_hit_info = &_list_of_mc_hit_info.at(loc);

        mc_hit_info->fMc   = mc;
        mc_hit_info->fFlag = chf;

        mc_hit_info->fType = mc->Type();

        int flagged_as_delta = chf->hasAnyProperty(StrawHitFlag::bkg);

        if (flagged_as_delta) {
          mc->fNShTaggedDelta += 1;
          fNHitsDeltaReco++;
        }

        if (ch_counted == 0) {
          mc->fListOfComboHits.push_back(ch);
          mc->fListOfChFlags.push_back(chf);

          if (flagged_as_delta) {
            mc->fNChTaggedDelta += 1;
          }
          ch_counted = 1;
        }

        loc++;
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

        if ((sim->pdgId() == PDGCode::e_minus) and (mc->NComboHits() >= _printElectronsMinNHits)) {

          float fr_sh = mc->NShTaggedDelta()/(mc->NSSCHits()+1.e-3);
          float fr_ch = mc->NChTaggedDelta()/(mc->NComboHits()+1.e-3);

          if (fr_ch < _printElectronsMaxFReco) {

            printf("* electron: sim.id: %6li mom: %8.3f time: %8.2f ",
                   sim->id().asInt(), mc->Momentum(), mc->Time());
            printf("nch:nchd:fr: %3i:%3i:%6.3f  nsh:nshd:fr %3i:%3i:%6.3f stations: %2i:%2i",
                   mc->NComboHits(),mc->NChTaggedDelta(),fr_ch,
                   mc->NSSCHits  (),mc->NShTaggedDelta(),fr_sh,
                   mc->fFirstStation, mc->fLastStation);
            printf(" \n");
//-----------------------------------------------------------------------------
// check if want to print hits - print 1-straw hits....
//-----------------------------------------------------------------------------
            if (_printElectronHits != 0) {
                                        // these are 'combo' ComboHits
              int nh  = mc->NComboHits();

              _hlp->printComboHit(0,0,"banner",0,0);

              for (int i=0; i<nh; i++) {
                const ComboHit* ch = mc->Ch(i);
                                        // flags here are the combo hit flags, not 1-straw hit flags
                int flags  = *((int*) mc->ChFlag(i));
                                        // use sim_id of the first SSCH ...
                int ind    = ch->indexArray()[0];

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
      printf("* ComboHits \n");
//-----------------------------------------------------------------------------
// print ComboHits
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_chCollTag.encode().data(),
                                    _chfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }

    if (_printSingleComboHits) {
      printf("* Single straw ComboHits tag:  %s\n",_shCollTag.encode().data());
//-----------------------------------------------------------------------------
// print ComboHits
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_shCollTag.encode().data(),
                                    _shfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }
  }

// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(DeltaFinderAna)

}
