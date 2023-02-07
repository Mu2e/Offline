#ifndef __CalPatRec_DeltaFinderDiag_hh__
#define __CalPatRec_DeltaFinderDiag_hh__

#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "Offline/CalPatRec/inc/HlPrint.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

using namespace std;

namespace mu2e {

  using namespace DeltaFinderTypes;

  // namespace DeltaFinderTypes {
  //   struct McPart_t;
  //   struct Data_t;
  // }

  // using DeltaFinderTypes::McPart_t;
  // using DeltaFinderTypes::Data_t;

  class SimParticle;
  class StrawDigiMCCollection;

  class DeltaFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets =  10,
      kNSeedHistSets  = 100,
      kNSeed2HistSets = 100,
      kNDeltaHistSets = 100,
      kNMcHistSets    = 200
    };

    struct SeedHist_t {
      TH1F*  fChi2TotN;
      TH1F*  fChi2ParN;
      TH1F*  fChi2PerpN;
      TH1F*  fChi2Delta;
      TH1F*  fChi2DeltaPar;
      TH1F*  fChi2DeltaPerp;
      TH1F*  fHitChi2Min;                        // chi2 of the first two hits along the wire
      TH1F*  fChi2Neighbour;
      TH1F*  fChi2Radial;
      // TH1F*  fNFacesWithHits;
      TH1F*  fNHitsPerFace;
      TH1F*  fNHitsPerSeed;
      TH1F*  fSeedRadius;
      TH1F*  fSeedMomentum;
      TH1F*  fNHits;
      TH1F*  fNSim;
      TH1F*  fNMom;
      TH1F*  fEDep;
      TH1F*  fDt;
    };

    struct Seed2Hist_t {
      TH1F*  fDt;                       // dt between the two seeds in teh same station
    };

    struct DeltaHist_t {
      TH1F*  fNHits;
      TH1F*  fNHitsCE;
      TH1F*  fNStrawHits;
      TH1F*  fNSeeds;
      TH1F*  fMcMom;
      TH1F*  fPDGCode;
      TH1F*  fDxy;
      TH1F*  fEDep;
      TH2F*  fNBestVsN;
      TH1F*  fFBest;
    };

    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom;
      TH1F*  fNHits;
      TH1F*  fNHitsDelta;
      TH1F*  fFractReco;
      TH1F*  fMaxSeg;
      TH2F*  fFractRecoVsNHits;
      TH1F*  fHitDt;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      //      TH1F*  fNSecondHits;
      TH1F*  fNSeeds;
      TH2F*  fNSeedsVsStation;
      TH1F*  fNMc;
      TH1F*  fPDGCode;
      TH1F*  fNMcHits;
      TH1F*  fNDelta;
      TH1F*  fNDeltaHitsT;
      TH1F*  fNDeltaHitsR;
    };

    struct Hist_t {
      EventHist_t* fEvent[kNEventHistSets];
      SeedHist_t*  fSeed [kNSeedHistSets ];
      Seed2Hist_t* fSeed2[kNSeed2HistSets];
      DeltaHist_t* fDelta[kNDeltaHistSets];
      McHist_t*    fMc   [kNMcHistSets   ];
    };
//-----------------------------------------------------------------------------
// additional parameters of the DeltaSeed , for analysis
//-----------------------------------------------------------------------------
    struct DeltaSeedPar_t {
      float dt;                           // time residual to delta
      McPart_t*    fMcPart[kNFaces];
      int          fNSim;                 // N different simID's cof the hits
      int          fNMom;                 // N(mothers) of the particles which produced the hits
      int          fNHitsCE;              // number of associated CE hits
      McPart_t*    fPreSeedMcPart[2];     // McPart_t's for initial stereo intersection

      int MCTruth () { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
    };

    std::vector<DeltaSeedPar_t>           _seedPar[kNStations];

  protected:

    bool                                  _mcDiag;
    int                                   _printOTracker;
    int                                   _printComboHits;
    int                                   _printElectrons;
    int                                   _printElectronsHits;
    int                                   _printElectronsMinNHits;
    float                                 _printElectronsMaxFReco;
    float                                 _printElectronsMinMom;
    float                                 _printElectronsMaxMom;
    int                                   _printDeltaSeeds;
    int                                   _printDeltaCandidates;
    int                                   _printShcol;
    int                                   _printSeedNParents;   // n(hits) for seeds with > 2 parents

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    int                                   _eventNumber;
    int                                   _nDeltaHitsTot;
    int                                   _nDeltaHitsReco;

    TObjArray                             _list_of_mc_particles; // list_of_particles with hits in the tracker
    TObjArray                             _list_of_mc_part_hit ; // for each StrawHit, pointer to its McPart

    Data_t*                               _data;                 // shared data, passed from the caller

    Hist_t                                _hist;

  public:

    DeltaFinderDiag(const fhicl::ParameterSet& PSet);
    DeltaFinderDiag(const fhicl::Table<mu2e::DeltaFinderTypes::Config>& config);
    ~DeltaFinderDiag();

  private:

    void        bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir);
    void        bookSeedHistograms (SeedHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookSeed2Histograms(Seed2Hist_t* Hist, art::TFileDirectory* Dir);
    void        bookDeltaHistograms(DeltaHist_t* Hist, art::TFileDirectory* Dir);
    void        bookMcHistograms   (McHist_t*    Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms(EventHist_t* Hist);
    void        fillSeedHistograms (SeedHist_t*  Hist, DeltaSeed*      Seed , DeltaSeedPar_t* SeedPar);
    void        fillSeed2Histograms(Seed2Hist_t* Hist, DeltaSeed*      Seed1, DeltaSeed* Seed2  );
    void        fillDeltaHistograms(DeltaHist_t* Hist, DeltaCandidate* Delta);
    void        fillMcHistograms   (McHist_t*    Hist, McPart_t*       Mc   );

    McPart_t*   findParticle (const SimParticle* Sim);
    int         InitMcDiag      ();
    int         associateMcTruth();

    void        printDeltaSeed(DeltaSeed* Seed, DeltaSeedPar_t* Sp, int PrintBanner);

    void        printHitData(const HitData_t* Hd, int PrintBanner);
    void        printOTracker();

    void        printComboHit(const ComboHit* Sh, int PrintBanner);
    void        printComboHitCollection();
//-----------------------------------------------------------------------------
// overriden virtual functions of the base class
//-----------------------------------------------------------------------------
  public:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1 ) override ;
    virtual int debug         (void* Data, int Mode = -1 ) override ;
  };


//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  DeltaFinderDiag::DeltaFinderDiag(const fhicl::ParameterSet& PSet):
    _mcDiag                (PSet.get<bool>         ("mcDiag"                       )),
    _printOTracker         (PSet.get<int>          ("printOTracker"                )),
    _printElectrons        (PSet.get<int>          ("printElectrons"               )),
    _printElectronsHits    (PSet.get<int>          ("printElectronsHits"           )),
    _printElectronsMinNHits(PSet.get<int>          ("printElectronsMinNHits"       )),
    _printElectronsMaxFReco(PSet.get<float>        ("printElectronsMaxFReco"       )),
    _printElectronsMinMom  (PSet.get<float>        ("printElectronsMinMom"         )),
    _printElectronsMaxMom  (PSet.get<float>        ("printElectronsMaxMom"         )),
    _printDeltaSeeds       (PSet.get<int>          ("printDeltaSeeds"              )),
    _printDeltaCandidates  (PSet.get<int>          ("printDeltaCandidates"         )),
    _printShcol            (PSet.get<int>          ("printShcol"                   )),
    _printSeedNParents     (PSet.get<int>          ("printSeedNParents"            ))
  {
    printf(" DeltaFinderDiag::DeltaFinderDiag : HOORAY PSet! \n");
    //    _firstCall   =  1;
    //    _timeOffsets = NULL;

    if (_mcDiag != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  DeltaFinderDiag::DeltaFinderDiag(const fhicl::Table<mu2e::DeltaFinderTypes::Config>& config):
    _mcDiag                (config().mcDiag()                ),
    _printOTracker         (config().printOTracker()         ),
    _printComboHits        (config().printComboHits()        ),
    _printElectrons        (config().printElectrons()        ),
    _printElectronsHits    (config().printElectronsHits()    ),
    _printElectronsMinNHits(config().printElectronsMinNHits()),
    _printElectronsMaxFReco(config().printElectronsMaxFReco()),
    _printElectronsMinMom  (config().printElectronsMinMom()  ),
    _printElectronsMaxMom  (config().printElectronsMaxMom()  ),
    _printDeltaSeeds       (config().printDeltaSeeds()       ),
    _printDeltaCandidates  (config().printDeltaCandidates()  ),
    _printShcol            (config().printShcol()            ),
    _printSeedNParents     (config().printSeedNParents()     )
  {
    printf(" DeltaFinderDiag::DeltaFinderDiag : HOORAY Config! \n");
    //    _firstCall   =  1;

    if (_mcDiag != 0) _mcUtils = art::make_tool  <McUtilsToolBase>(config().mcUtils.get_PSet());
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  DeltaFinderDiag::~DeltaFinderDiag() {
  }

//-----------------------------------------------------------------------------
  McPart_t* DeltaFinderDiag::findParticle(const SimParticle* Sim) {
    McPart_t* found(0);

    int n = _list_of_mc_particles.GetEntriesFast();

    for (int i=0; i<n; i++) {
      McPart_t* mc = (McPart_t*) _list_of_mc_particles.At(i);
      if (mc->fSim == Sim) {
        found = mc;
        break;
      }
    }

    return found;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookDeltaHistograms(DeltaHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fNSeeds     = Dir->make<TH1F>("nseeds"    , "N(seeds)"        ,  20, 0.,  20.);
    Hist->fNHits      = Dir->make<TH1F>("nch"       , "N(combo hits)"   , 200, 0., 200.);
    Hist->fNHitsCE    = Dir->make<TH1F>("nch_ce"    , "N(combo hits CE)",  50, 0.,  50.);
    Hist->fNStrawHits = Dir->make<TH1F>("nsh"       , "N(straw hits)"   , 200, 0., 200.);
    Hist->fPDGCode    = Dir->make<TH1F>("pdg"       , "PDG code"        ,1000, -2500., 2500.);
    Hist->fMcMom      = Dir->make<TH1F>("mc_mom"    , "N(hits)"         , 200, 0., 200.);
    Hist->fDxy        = Dir->make<TH1F>("dxy"       , "Delta Dxy"       , 500, 0., 200.);
    Hist->fEDep       = Dir->make<TH1F>("edep"      , "mean E(Dep)"     , 100, 0., 0.01);
    Hist->fNBestVsN   = Dir->make<TH2F>("nbest_vs_n", "N(bset) vs N"    , 100, 0., 100.,100,0,100);
    Hist->fFBest      = Dir->make<TH1F>("fbest"     , "N(best)/N"       , 100, 0.,    1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number", 100, 0., 100000.);
    //    Hist->fNSecondHits     = Dir->make<TH1F>("nhit2" , "N(second hits)", 100, 0., 100.);
    Hist->fNSeeds          = Dir->make<TH1F>("nseeds", "N(seeds)"   , 200, 0., 2000.);
    Hist->fNSeedsVsStation = Dir->make<TH2F>("ns_vs_st", "N(seeds) vs station", 20, 0., 20.,200,0,200);

    Hist->fNMc             = Dir->make<TH1F>("nmc"       , "N(MC particles)", 100, 0., 1000.);
    Hist->fPDGCode         = Dir->make<TH1F>("pdg_code"  , "PDG Code"       ,1000, -2500., 2500.);
    Hist->fNMcHits         = Dir->make<TH1F>("n_mc_hits" , "N(MC hits)"     , 250, 0., 500.);
    Hist->fNDelta          = Dir->make<TH1F>("n_delta"   , "N(reco deltas)" , 250, 0., 250.);
    Hist->fNDeltaHitsT     = Dir->make<TH1F>("n_delta_ht", "N(delta hits T)", 500, 0., 5000.);
    Hist->fNDeltaHitsR     = Dir->make<TH1F>("n_delta_hr", "N(delta hits R)", 500, 0., 5000.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookMcHistograms(McHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fPDGCode    = Dir->make<TH1F>("pdg"  , "PDG code"        ,1000, -2500, 2500);
    Hist->fMom        = Dir->make<TH1F>("mom"  , "momentum"        , 200, 0., 200.);
    Hist->fNHits      = Dir->make<TH1F>("nhits", "N(hits)"         , 200, 0., 200.);
    Hist->fNHitsDelta = Dir->make<TH1F>("nhitsr", "N(hits reco)"   , 200, 0., 200.);
    Hist->fFractReco  = Dir->make<TH1F>("fractr", "NR/N"           , 100, 0.,   1.);
    Hist->fMaxSeg     = Dir->make<TH1F>("max_seg", "Max N Segments",  20, 0.,  20.);
    Hist->fHitDt      = Dir->make<TH1F>("hit_dt" , "Hit TMax-TMin" , 100, 0., 200.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookSeedHistograms(SeedHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fChi2TotN       = Dir->make<TH1F>("chi2totn"       , "Chi2Tot/N"         , 1000, 0., 100.);
    Hist->fChi2ParN       = Dir->make<TH1F>("chi2parn"       , "Chi2Par/N"         , 1000, 0., 100.);
    Hist->fChi2PerpN      = Dir->make<TH1F>("chi2perpn"      , "Chi2Perp/N"        , 1000, 0., 100.);
    Hist->fChi2Delta      = Dir->make<TH1F>("chi2_delta"     , "Chi2Delta"         , 1000, 0., 100.);
    Hist->fChi2DeltaPar   = Dir->make<TH1F>("chi2_delta_par" , "Chi2DeltaPar"      , 1000, 0., 100.);
    Hist->fChi2DeltaPerp  = Dir->make<TH1F>("chi2_delta_perp", "Chi2DeltaPerp"     , 1000, 0., 100.);
    Hist->fHitChi2Min     = Dir->make<TH1F>("hitchi2min"     , "Hit Chi (min)"     , 1000, 0.,  50.);
    Hist->fChi2Neighbour  = Dir->make<TH1F>("chi2_nb"        , "Chi2 neighbour"    , 1000, 0.,  50.);
    Hist->fChi2Radial     = Dir->make<TH1F>("chi2_r"         , "Chi2 radial"       , 1000, 0.,  50.);
    Hist->fNHitsPerFace   = Dir->make<TH1F>("nhits_face"     , "Number of hits per face", 20, 0., 20.);
    Hist->fNHits          = Dir->make<TH1F>("nhits"          , "Number of hits"      , 10, 0., 10.);
    Hist->fSeedRadius     = Dir->make<TH1F>("rad"            , "Seed radius", 500, 0., 1000.);
    Hist->fSeedMomentum   = Dir->make<TH1F>("mom"            , "Seed momentum", 400, 0., 400.);
    Hist->fNSim           = Dir->make<TH1F>("nsim"           , "N(sim particls)"     , 10,    0.,  10.);
    Hist->fNMom           = Dir->make<TH1F>("nmom"           , "N(MC moms)"          , 10,    0.,  10.);
    Hist->fEDep           = Dir->make<TH1F>("edep"           , "mean E(Dep)"         ,100,    0.,   0.01);
    Hist->fDt             = Dir->make<TH1F>("dt"             , "delta(T)"            ,200, -100., 100.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookSeed2Histograms(Seed2Hist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fDt              = Dir->make<TH1F>("dt"          , "delta(T)"            ,100, -100, 100);
  }

//-----------------------------------------------------------------------------
// this routine is called once per job (likely, from beginJob)
// TH1::AddDirectory makes sure one can have histograms with the same name
// in different subdirectories
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    TH1::AddDirectory(0);
    char folder_name[20];
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book seed histograms
//-----------------------------------------------------------------------------
    int book_seed_histset[kNSeedHistSets];
    for (int i=0; i<kNSeedHistSets; i++) book_seed_histset[i] = 0;

    book_seed_histset[ 0] = 1;          // all seeds
    book_seed_histset[ 1] = 1;          // nhits=1
    book_seed_histset[ 2] = 1;          // nhits=2
    book_seed_histset[ 3] = 1;          // nhits=3
    book_seed_histset[ 4] = 1;          // nhits=4

    book_seed_histset[11] = 1;          // MC-2-true, protons + deut
    book_seed_histset[12] = 1;          // MC-2-true, e-   0 < p <  20
    book_seed_histset[13] = 1;          // MC-2-true, e-  20 < p <  80
    book_seed_histset[14] = 1;          // MC-2-true, e-  80 < p < 110
    book_seed_histset[15] = 1;          // MC-2-true, e- 110 < p < inf
    book_seed_histset[16] = 1;          // MC-2-true, positrons
    book_seed_histset[17] = 1;          // MC-2-true, muons
    book_seed_histset[18] = 1;          // MC-2-true, pions

    book_seed_histset[20] = 1;          // all seeds from deltas

    book_seed_histset[21] = 1;          // delta(nhits>=5) seed nhits=1
    book_seed_histset[22] = 1;          // delta(nhits>=5) seed nhits=2
    book_seed_histset[23] = 1;          // delta(nhits>=5) seed nhits=3
    book_seed_histset[24] = 1;          // delta(nhits>=5) seed nhits=4

    book_seed_histset[31] = 1;          // delta(nseeds>3) seed nhits=1
    book_seed_histset[32] = 1;          // delta(nseeds>3) seed nhits=2
    book_seed_histset[33] = 1;          // delta(nseeds>3) seed nhits=3
    book_seed_histset[34] = 1;          // delta(nseeds>3) seed nhits=4

    book_seed_histset[90] = 1;          // truth = 0
    book_seed_histset[91] = 1;          // N(CE hits) > 0

    for (int i=0; i<kNSeedHistSets; i++) {
      if (book_seed_histset[i] != 0) {
        sprintf(folder_name,"seed_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fSeed[i] = new SeedHist_t;
        bookSeedHistograms(_hist.fSeed[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book seed2 histograms
//-----------------------------------------------------------------------------
    int book_seed2_histset[kNSeed2HistSets];
    for (int i=0; i<kNSeed2HistSets; i++) book_seed2_histset[i] = 0;

    book_seed2_histset[ 0] = 1;          // all pairs of N(hits)>=2 seeds
    book_seed2_histset[ 1] = 1;          // all pairs of seeds found in a search mode

    for (int i=0; i<kNSeed2HistSets; i++) {
      if (book_seed2_histset[i] != 0) {
        sprintf(folder_name,"seed2_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fSeed2[i] = new Seed2Hist_t;
        bookSeed2Histograms(_hist.fSeed2[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book delta histograms
//-----------------------------------------------------------------------------
    int book_delta_histset[kNDeltaHistSets];
    for (int i=0; i<kNDeltaHistSets; i++) book_delta_histset[i] = 0;

    book_delta_histset[ 0] = 1;                // all  deltas
    book_delta_histset[ 1] = 1;                // prot+deut
    book_delta_histset[ 2] = 1;                // low-mom e-
    book_delta_histset[ 6] = 1;                // e+
    book_delta_histset[ 7] = 1;                // muons
    book_delta_histset[ 8] = 1;                // pions

    book_delta_histset[10] = 1;                // N(hits)>=5
    book_delta_histset[11] = 1;                // N(hits)>=5 prot+deut
    book_delta_histset[12] = 1;                // N(hits)>=5 low-mom e-
    book_delta_histset[16] = 1;                // N(hits)>=5 e+
    book_delta_histset[17] = 1;                // N(hits)>=5 muons
    book_delta_histset[18] = 1;                // N(hits)>=5 pions

    book_delta_histset[20] = 1;                // N(hits)>=5 and NHitsCE>2

    for (int i=0; i<kNDeltaHistSets; i++) {
      if (book_delta_histset[i] != 0) {
        sprintf(folder_name,"delta_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fDelta[i] = new DeltaHist_t;
        bookDeltaHistograms(_hist.fDelta[i],&dir);
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
    book_mc_histset[  4] = 1;                // electrons fTime > 550 with last>first and 4+ hits
    book_mc_histset[  5] = 1;                // electrons fTime > 550 with last>first, 5+ hits, and reco delta
    book_mc_histset[  6] = 1;                // electrons fTime > 550 with last>first, 5+ hits, and p < 20

    book_mc_histset[100] = 1;                // electrons fTime > 550 and 20 < p < 80 MeV/c
    book_mc_histset[101] = 1;                // electrons fTime > 550 and p > 80 MeV/c

    for (int i=0; i<kNMcHistSets; i++) {
      if (book_mc_histset[i] != 0) {
        sprintf(folder_name,"mc_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fMc[i] = new McHist_t;
        bookMcHistograms(_hist.fMc[i],&dir);
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillDeltaHistograms(DeltaHist_t* Hist, DeltaCandidate* Delta) {
    int n_seeds = Delta->fNSeeds;
    Hist->fNSeeds->Fill(n_seeds);
    Hist->fNHits->Fill(Delta->fNHits);
    Hist->fNHitsCE->Fill(Delta->fNHitsCE);
    Hist->fNStrawHits->Fill(Delta->fNStrawHits);

    float mom(-1), pdg_code(-1.e6);
    if (Delta->fMcPart) {
      mom      = Delta->fMcPart->Momentum();
      pdg_code = Delta->fMcPart->fPdgID;
    }

    Hist->fPDGCode->Fill(pdg_code);
    Hist->fMcMom->Fill(mom);
    Hist->fEDep->Fill(Delta->EDep());
    Hist->fNBestVsN->Fill(Delta->NHits(),Delta->NHitsMcP());

    float fbest = float(Delta->NHitsMcP())/Delta->NHits();
    Hist->fFBest->Fill(fbest);

    for(int is=0; is<kNStations; ++is) {
      Hist->fDxy->Fill(Delta->dxy[is]);
    }
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();

    Hist->fEventNumber->Fill(event_number);
    Hist->fNSeeds->Fill(_data->nSeedsTot());

    for (int is=0; is<kNStations; ++is) {

      int nseeds = _data->NSeeds(is);
      Hist->fNSeedsVsStation->Fill(is,nseeds);

      // for (int f=0; f<kNFaces-1; ++f) { // loop over faces
      //   for (int p=0; p<kNPanelsPerFace; ++p) { // loop over panels
      //     PanelZ_t* panelz = &_data->oTracker[is][f][p];

      //     // int nhits =  panelz->fHitData.size();
      //     // for (int i=0; i<nhits; i++) {
      //     //   int counter  = panelz->fHitData.at(i).fNSecondHits;
      //     //   Hist->fNSecondHits->Fill(counter);
      //     // }
      //   }
      // }
    }

    int ndelta = _data->nDeltaCandidates();
    Hist->fNDelta->Fill(ndelta);
//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.GetEntriesFast();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _list_of_mc_particles.At(i);
      int n_mc_hits = mc->fListOfHits.size();

      Hist->fNMcHits->Fill(n_mc_hits);
      Hist->fPDGCode->Fill(mc->fPdgID);
    }

    Hist->fNMc->Fill(nmc);
    Hist->fNDeltaHitsT->Fill(_nDeltaHitsTot);
    Hist->fNDeltaHitsR->Fill(_nDeltaHitsReco);
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillSeedHistograms(SeedHist_t* Hist, DeltaSeed* Seed, DeltaSeedPar_t* Par) {

    Intersection_t res;

    Hist->fChi2TotN->Fill(Seed->Chi2TotN());
    Hist->fChi2PerpN->Fill(Seed->Chi2PerpN());
    Hist->fChi2ParN ->Fill(Seed->Chi2ParN ());
    Hist->fChi2Delta->Fill(Seed->Chi2Delta());
    Hist->fChi2DeltaPar->Fill(Seed->Chi2DeltaPar());
    Hist->fChi2DeltaPerp->Fill(Seed->Chi2DeltaPerp());

    int face0 = Seed->SFace(0);
    int face1 = Seed->SFace(1);

    if (face1 >= 0) {
//-----------------------------------------------------------------------------
// full-fledged stereo seed
//-----------------------------------------------------------------------------
      const HitData_t* hd1 = Seed->HitData(face0);
      const HitData_t* hd2 = Seed->HitData(face1);

      DeltaFinderTypes::findIntersection(hd1,hd2,&res);

      double chi2_1 = res.wd1*res.wd1/hd1->fSigW2;
      double chi2_2 = res.wd2*res.wd2/hd2->fSigW2;

      Hist->fHitChi2Min->Fill(chi2_1);
      Hist->fHitChi2Min->Fill(chi2_2);

      // int nh1 = hd1->fHit->nStrawHits();
      // int nh2 = hd2->fHit->nStrawHits();

      double _sigmaR = 10;

      for (int face=0; face<kNFaces; face++) {
//-----------------------------------------------------------------------------
// skip seed faces
//-----------------------------------------------------------------------------
        if ((face == face0) || (face == face1))                       continue;

        const HitData_t* hd  = Seed->HitData(face);
        if (hd == nullptr)                                            continue;
//-----------------------------------------------------------------------------
// reproduce DeltaFinder algorithm
// *FIXME* - replace with the call to the algorithm function...
//-----------------------------------------------------------------------------
        XYZVectorF dxyz = hd->fHit->pos()-Seed->CofM;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
        const XYZVectorF& wdir     = hd->fHit->wdir();
        XYZVectorF        d_par    = (dxyz.Dot(wdir))/(wdir.Dot(wdir))*wdir;
        XYZVectorF        d_perp_z = dxyz-d_par;
        float  d_perp              = d_perp_z.rho();
        // double sigw                = hd->fSigW;
        float  chi2_par            = (d_par.R()*d_par.R())/hd->fSigW2;
        float  chi2_perp           = (d_perp/_sigmaR)*(d_perp/_sigmaR);
        float  chi2r               = chi2_par + chi2_perp;
        Hist->fChi2Radial->Fill(chi2r);
      }
    }

    // Hist->fNFacesWithHits->Fill(Seed->fNFacesWithHits);
    Hist->fSeedRadius->Fill    (Seed->CofM.rho());

    double mom (-1.);
    if (Par->fPreSeedMcPart[0]) mom = Par->fPreSeedMcPart[0]->Momentum();

    Hist->fSeedMomentum->Fill(mom);

    Hist->fNHits->Fill(Seed->NHits());
    Hist->fNSim->Fill(Par->fNSim);
    Hist->fNMom->Fill(Par->fNMom);

    Hist->fEDep->Fill(Seed->EDep());

    Hist->fDt->Fill(Par->dt);
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillSeed2Histograms(Seed2Hist_t* Hist, DeltaSeed* Seed1, DeltaSeed* Seed2) {
    float dt = Seed1->TMean()-Seed2->TMean();
    Hist->fDt->Fill(dt);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillMcHistograms(McHist_t* Hist, McPart_t* Mc) {

    float mom = Mc->Momentum();

    Hist->fPDGCode->Fill(Mc->fPdgID);
    Hist->fMom->Fill(mom);
    Hist->fNHits->Fill(Mc->NHits());
    Hist->fNHitsDelta->Fill(Mc->fNHitsDelta);

    float freco = Mc->fNHitsDelta/(Mc->NHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNHits->Fill(Mc->NHits(),freco);

    int max_nseg = Mc->fLastStation-Mc->fFirstStation+1;
    Hist->fMaxSeg->Fill(max_nseg);

    float dt = Mc->HitDt();
    Hist->fHitDt->Fill(dt);
  }

//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// start from precalculating MC-specific info
//-----------------------------------------------------------------------------
    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
        _eventNumber = en;
        InitMcDiag();
        associateMcTruth();
      }
    }
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// fill seed histograms
//-----------------------------------------------------------------------------
    for (int s=0; s<kNStations; ++s) {
      int nseeds = _data->NSeeds(s);
      for(int se=0; se<nseeds; ++se) {
        DeltaSeed*      seed = _data->deltaSeed(s,se);
        DeltaSeedPar_t* sp   = &_seedPar[s].at(se);

        float mom    = sp->fPreSeedMcPart[0]->Momentum();
        int   pdg_id = sp->fPreSeedMcPart[0]->fPdgID;

        sp->dt = -1.e6;
        if (seed->fDeltaIndex >= 0) {
          DeltaCandidate* dc = _data->deltaCandidate(seed->fDeltaIndex);
          sp->dt = seed->TMean() - dc->T0(s);
        }

        fillSeedHistograms(_hist.fSeed[0],seed,sp);

        if      (seed->NHits() == 1) fillSeedHistograms(_hist.fSeed[ 1],seed,sp);
        else if (seed->NHits() == 2) fillSeedHistograms(_hist.fSeed[ 2],seed,sp);
        else if (seed->NHits() == 3) fillSeedHistograms(_hist.fSeed[ 3],seed,sp);
        else if (seed->NHits() == 4) fillSeedHistograms(_hist.fSeed[ 4],seed,sp);

        if (sp->MCTruth()) {
//-----------------------------------------------------------------------------
// a delta seed which started by a stereo of two hits produced by the same particle
//-----------------------------------------------------------------------------
          if      (pdg_id >      2000) fillSeedHistograms(_hist.fSeed[11],seed,sp);
          else if (pdg_id ==  11) {
            if      (mom <   20)       fillSeedHistograms(_hist.fSeed[12],seed,sp);
            else if (mom <   80)       fillSeedHistograms(_hist.fSeed[13],seed,sp);
            else if (mom <  110)       fillSeedHistograms(_hist.fSeed[14],seed,sp);
            else                       fillSeedHistograms(_hist.fSeed[15],seed,sp);
          }
          else if (pdg_id      == -11) fillSeedHistograms(_hist.fSeed[16],seed,sp);
          else if (abs(pdg_id) ==  13) fillSeedHistograms(_hist.fSeed[17],seed,sp);
          else if (abs(pdg_id) == 211) fillSeedHistograms(_hist.fSeed[18],seed,sp);
        }
        else {
//-----------------------------------------------------------------------------
// perhaps fake seed - made out of hits produced by two different particles
//-----------------------------------------------------------------------------
          fillSeedHistograms(_hist.fSeed[90],seed,sp);
        }

        if (sp->fNHitsCE > 0) {
          fillSeedHistograms(_hist.fSeed[91],seed,sp);
          if (seed->Chi2TotN() > 50) {
            printf("* ERROR seed with chi2TotN  > 50\n");
            _data->printDeltaSeed(seed);
          }
        }

        DeltaCandidate* dc(nullptr);
        if (seed->fDeltaIndex >= 0) dc = _data->deltaCandidate(seed->fDeltaIndex);
        if (dc) {
//-----------------------------------------------------------------------------
// this seed has been associated with a delta
//-----------------------------------------------------------------------------
          if (dc->NHits() >= 5) {
            if      (seed->NHits() == 1) fillSeedHistograms(_hist.fSeed[21],seed,sp);
            else if (seed->NHits() == 2) fillSeedHistograms(_hist.fSeed[22],seed,sp);
            else if (seed->NHits() == 3) fillSeedHistograms(_hist.fSeed[23],seed,sp);
            else if (seed->NHits() == 4) fillSeedHistograms(_hist.fSeed[24],seed,sp);
          }

          if (dc->NSeeds() > 3) {
            if      (seed->NHits() == 1) fillSeedHistograms(_hist.fSeed[31],seed,sp);
            else if (seed->NHits() == 2) fillSeedHistograms(_hist.fSeed[32],seed,sp);
            else if (seed->NHits() == 3) fillSeedHistograms(_hist.fSeed[33],seed,sp);
            else if (seed->NHits() == 4) fillSeedHistograms(_hist.fSeed[34],seed,sp);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// fill seed2 histograms
// [0] : all pairs of nhits>=2 seeds
// [0] : all pairs of seeds found in a search mode  (the ones which potentially need pruning)
//-----------------------------------------------------------------------------
    for (int s=0; s<kNStations; ++s) {
      int nseeds = _data->NSeeds(s);

      for (int i1=0; i1<nseeds-1; i1++) {
        DeltaSeed* seed1 = _data->deltaSeed(s,i1);
        if (seed1->NHits() < 2) continue;

        for (int i2=i1+1; i2<nseeds; i2++) {

          DeltaSeed* seed2 = _data->deltaSeed(s,i2);
          if (seed2->NHits() < 2) continue;

          fillSeed2Histograms(_hist.fSeed2[0],seed1,seed2);
        }
      }

      for(int i1=0; i1<nseeds-1; i1++) {
        DeltaSeed* seed1 = _data->deltaSeed(s,i1);
        if (seed1->SFace(1) < 0) continue;

        for (int i2=i1+1; i2<nseeds; i2++) {

          DeltaSeed* seed2 = _data->deltaSeed(s,i2);
          if (seed2->SFace(1) < 0) continue;

          fillSeed2Histograms(_hist.fSeed2[1],seed1,seed2);
        }
      }
    }
//-----------------------------------------------------------------------------
// fill delta histograms
//-----------------------------------------------------------------------------
    int ndelta = _data->nDeltaCandidates();
    for (int i=0; i<ndelta; i++) {
      DeltaCandidate* delta = _data->deltaCandidate(i);
      int pdg_id = delta->fMcPart->fPdgID;
      float mom  = delta->fMcPart->Momentum();

      fillDeltaHistograms(_hist.fDelta[0],delta);

      if      (pdg_id > 2000)    fillDeltaHistograms(_hist.fDelta[1],delta);
      else if (pdg_id ==  11) {
        if (mom < 20)            fillDeltaHistograms(_hist.fDelta[2],delta);
      }
      else if (pdg_id      == -11)    fillDeltaHistograms(_hist.fDelta[6],delta);
      else if (abs(pdg_id) ==  13)    fillDeltaHistograms(_hist.fDelta[7],delta);
      else if (abs(pdg_id) == 211)    fillDeltaHistograms(_hist.fDelta[8],delta);

      if (delta->NHits() >= 5) {
        fillDeltaHistograms(_hist.fDelta[10],delta);

        if      (pdg_id      > 2000) fillDeltaHistograms(_hist.fDelta[11],delta);
        else if (pdg_id ==  11) {
          if (mom < 20)              fillDeltaHistograms(_hist.fDelta[12],delta);
        }
        else if (pdg_id      == -11) fillDeltaHistograms(_hist.fDelta[16],delta);
        else if (abs(pdg_id) ==  13) fillDeltaHistograms(_hist.fDelta[17],delta);
        else if (abs(pdg_id) == 211) fillDeltaHistograms(_hist.fDelta[18],delta);

        if (delta->fNHitsCE > 2)     fillDeltaHistograms(_hist.fDelta[20],delta);
      }
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.GetEntriesFast();

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _list_of_mc_particles.UncheckedAt(i);

      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if (mc->fPdgID == PDGCode::e_minus) {
        fillMcHistograms(_hist.fMc[1],mc);

        if (mc->Time() > 550) {
          fillMcHistograms(_hist.fMc[2],mc);

          if (mc->fLastStation > mc->fFirstStation) {
            fillMcHistograms(_hist.fMc[3],mc);

            if (mc->NHits() >= 5) {
              fillMcHistograms(_hist.fMc[4],mc);

              if (mc->fDelta != NULL) {
                fillMcHistograms(_hist.fMc[5],mc);
              }

              if (mc->Momentum() < 20) {
                fillMcHistograms(_hist.fMc[6],mc);
              }
            }
          }

          if ((mc->Momentum() > 20) && (mc->Momentum() < 80)) {
            fillMcHistograms(_hist.fMc[100],mc);
          }

          if (mc->Momentum() > 80) {
            fillMcHistograms(_hist.fMc[101],mc);
          }
        }
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker
// the hope is that it is shorter than the list of all particles
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::InitMcDiag() {
//-----------------------------------------------------------------------------
// memory cleanup after previous event
// assume that MC collections have been initialized
//-----------------------------------------------------------------------------
    int n = _list_of_mc_particles.GetEntriesFast();
    for (int i=0; i<n; i++) {
      McPart_t* p = (McPart_t*) _list_of_mc_particles.UncheckedAt(i);
      delete p;
    }

    _list_of_mc_particles.Clear();
    _list_of_mc_part_hit.Clear();

    int   nch           = _data->chcol->size();
    if (nch <= 0) return -1;

    const ComboHit* ch0 = &_data->chcol->at(0);

    _list_of_mc_part_hit.Expand(nch);

    for (int ist=0; ist<kNStations; ist++) {
      for (int face=0; face<kNFaces; face++) {
        FaceZ_t* fz = &_data->fFaceData[ist][face];

        int nhits = fz->fHitData.size();
        for (int ih=0; ih<nhits; ih++) {
          HitData_t*      hd = &fz->fHitData[ih];
          const ComboHit* ch = hd->fHit;
          size_t ich         = ch-ch0;  // hit index in the collection
//-----------------------------------------------------------------------------
// get MC particle associated with the particle hitting the first straw
// a combo hit could be made out of the straw hits produced by different particles
// in a neighbor straws - ignore that for now
//-----------------------------------------------------------------------------
          int loc = ch->index(0);
          const mu2e::SimParticle* sim = _mcUtils->getSimParticle(_data->event,loc);
//-----------------------------------------------------------------------------
// search if this particle has already been registered
//-----------------------------------------------------------------------------
          McPart_t* mc = findParticle(sim);

          if (mc == NULL) {
                                        // add new particle
            mc = new McPart_t(sim);
            _list_of_mc_particles.Add(mc);
            mc->fID       = _mcUtils->getID(sim);
            mc->fMotherID = _mcUtils->getMotherID(sim);
            mc->fPdgID    = _mcUtils->getPdgID(sim);
            mc->fStartMom = _mcUtils->getStartMom(sim);
          }
//-----------------------------------------------------------------------------
// list of hits produced by the particle
//-----------------------------------------------------------------------------
          mc->fListOfHits.push_back(hd);

          int station        = ch->strawId().station();

          if (station < mc->fFirstStation) mc->fFirstStation = station;
          if (station > mc->fLastStation ) mc->fLastStation  = station;

          if (ch->correctedTime() < mc->fTime   ) mc->fTime    = ch->correctedTime();
          if (ch->correctedTime() < mc->fHitTMin) mc->fHitTMin = ch->correctedTime();
          if (ch->correctedTime() > mc->fHitTMax) mc->fHitTMax = ch->correctedTime();
//-----------------------------------------------------------------------------
// list of MC particles parallel to StrawHitCollection
//-----------------------------------------------------------------------------
          _list_of_mc_part_hit.AddAt(mc,ich);
        }
      }
    }

    if (_data->debugLevel > 10) {
      int nmc = _list_of_mc_particles.GetEntriesFast();
      printf(" N(MC particles with hits in the tracker: %5i\n",nmc);
      printf("    i     SimID        PdgID  NHits   Momentum  Time   FirstSt LastSt\n");
      for (int i=0; i<nmc; i++) {
        McPart_t* mc = (McPart_t*)_list_of_mc_particles.At(i);
        if (mc) {
          printf(" %4i  %10i %10i  %5i %10.3f %8.1f %5i %5i\n",
                 i,mc->fID,mc->fPdgID,
                 mc->NHits(),
                 mc->Momentum(),
                 mc->Time(),
                 mc->fFirstStation,mc->fLastStation);
        }
        else {
          printf(" %4i  mc = nullptr\n",i);
        }
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
// for each DeltaSeed, create a list of SimParticle*'s parallel to its list of combo hits
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::associateMcTruth() {

    if (_data->chcol->size() == 0) return 0;

    const ComboHit* hit0 = &_data->chcol->at(0);

    for (int is=0; is<kNStations; is++) {
      int nseeds = _data->NSeeds(is);

      _seedPar[is].resize(nseeds);

      for (int i=0; i<nseeds; ++i) {
        DeltaSeed* seed = _data->deltaSeed(is,i);

        DeltaSeedPar_t* sp = &_seedPar[is].at(i);
//-----------------------------------------------------------------------------
// define MC pointers (SimParticle's) for the first two, "pre-seed", hits
//-----------------------------------------------------------------------------
        sp->fNHitsCE          = 0;
        sp->fPreSeedMcPart[0] = NULL;
        sp->fPreSeedMcPart[1] = NULL;

        int face0 = seed->SFace(0);
        if (seed->HitData(face0) != nullptr) {
          const HitData_t* hd1 = seed->HitData(face0);
          int loc1 = hd1->fHit-hit0;
          sp->fPreSeedMcPart[0] = (McPart_t*) _list_of_mc_part_hit.At(loc1);
        }

        int face1 = seed->SFace(1);
        if (face1 >= 0) {
          const HitData_t* hd2 = seed->HitData(face1);
          int loc2 = hd2->fHit-hit0;
          sp->fPreSeedMcPart[1] = (McPart_t*) _list_of_mc_part_hit.At(loc2);
        }
//-----------------------------------------------------------------------------
// define MC pointers (McPart_t's) for hits in all faces
// assume one CombiHit per face
//-----------------------------------------------------------------------------
        for (int face=0; face<kNFaces; face++) {
          const HitData_t* hd = seed->HitData(face);
          sp->fMcPart[face] = nullptr;
          if (hd == nullptr)                                          continue;
          int loc           = hd->fHit-hit0;         // ComboHit index in the collection
          McPart_t* mc      = (McPart_t*) _list_of_mc_part_hit.At(loc);
          sp->fMcPart[face] = mc;
//-----------------------------------------------------------------------------
// count CE hits
//-----------------------------------------------------------------------------
          if ((mc->fPdgID == PDGCode::e_minus) && (mc->fStartMom > 95) && (mc->fStartMom <110)) {
            sp->fNHitsCE += 1;
            if (seed->fDeltaIndex >= 0) {
              DeltaCandidate* dc = _data->deltaCandidate(seed->fDeltaIndex);
              dc->fNHitsCE += 1;
            }
          }
        }
//-----------------------------------------------------------------------------
// count number of different simID's contributing
//-----------------------------------------------------------------------------
        sp->fNSim  = 0;
        sp->fNMom = 0;
        int simid[4], nsim(0), momid[4], nmom(0);
        for (int face=0; face<kNFaces; face++) {
          McPart_t* mc = sp->fMcPart[face];
          if (mc == nullptr) continue;

          int found = 0;
          for (int i=0; i<nsim; i++) {
            if (mc->fID == simid[i]) {
              found = 1;
              break;
            }
          }

          if (found == 0) {
            simid[nsim] = mc->fID;
            nsim += 1;
          }
//-----------------------------------------------------------------------------
// repeat with mothers
//-----------------------------------------------------------------------------
          found = 0;
          for (int i=0; i<nmom; i++) {
            if (mc->fMotherID == momid[i]) {
              found = 1;
              break;
            }
          }

          if (found == 0) {
            momid[nmom] = mc->fMotherID;
            nmom += 1;
          }
        }

        sp->fNSim = nsim;
        sp->fNMom = nmom;
//-----------------------------------------------------------------------------
// print potentially problematic cases
//-----------------------------------------------------------------------------
        if (_printSeedNParents > 0) {
          if ((seed->fNHits == _printSeedNParents) and (sp->fNSim > 1)) {
            printf("  >>> WARNING: suspisious DeltaSeed: %i hits, Nmc = %i\n",seed->fNHits,sp->fNSim);
            printDeltaSeed(seed,sp,0);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// proceed with the reconstructed delta candidates
//-----------------------------------------------------------------------------
    int const    max_part(1000);
    McPart_t*    part [max_part];
    int          nhits[max_part];

    int ndelta = _data->nDeltaCandidates();
    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = _data->deltaCandidate(idelta);
      int npart          = 0;

      for (int is=dc->fFirstStation; is<=dc->fLastStation; is++) {
        DeltaSeed* ds = dc->Seed(is);
        if (ds == 0) continue;

        DeltaSeedPar_t* sp = &_seedPar[is].at(ds->Index());

        for (int face=0; face<kNFaces; face++) {
//-----------------------------------------------------------------------------
// assign delta index to each hit flagged as delta
//-----------------------------------------------------------------------------
          HitData_t* hd = (HitData_t*) ds->HitData(face);
          if (hd == nullptr)                                          continue;
          hd->fDeltaIndex = idelta;
//-----------------------------------------------------------------------------
// try to identify a reconstructed delta-electron with the MC particle
//-----------------------------------------------------------------------------
          McPart_t* mc = sp->fMcPart[face];  // list parallel to the list of hits

          int found(0);

          for (int ip=0; ip<npart; ip++) {
            if (mc == part[ip]) {
              found = 1;
              nhits[ip]++;
              break;
            }
          }

          if (found == 0) {
            if (npart < max_part) {
              part [npart] = mc;
              nhits[npart] = 1;
              npart += 1;
            }
            else {
              printf("DeltaFinder::associateMcTruth ERROR: npart >= max_part (%i)\n",max_part);
            }
          }
        }
      }
//-----------------------------------------------------------------------------
// look at the particles contributed to the seed and determine
// the "best" one - the one contributed the most number of hits
//-----------------------------------------------------------------------------
      int max_hits(-1), ibest(-1);

      for (int i=0; i<npart; i++) {
        if (nhits[i] > max_hits) {
          max_hits = nhits[i];
          ibest    = i;
        }
      }

      if (ibest >= 0) {
        dc->fMcPart         = part [ibest];
        dc->fNHitsMcP       = nhits[ibest];
        dc->fMcPart->fDelta = dc;
      }
    }
//-----------------------------------------------------------------------------
// now, for each MC electron calculate the number of hits flagged as delta
//-----------------------------------------------------------------------------
    int nmc    = _list_of_mc_particles.GetEntriesFast();

    _nDeltaHitsTot  = 0;
    _nDeltaHitsReco = 0;

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _list_of_mc_particles.At(i);
      mc->fNHitsDelta = 0;
//-----------------------------------------------------------------------------
// loop over the hits of MC delta electron and calculate fraction of them which have
// been reconstructed as hits of reconstructed delta electrons
//-----------------------------------------------------------------------------
      int nh = mc->fListOfHits.size();
      for (int ih=0; ih<nh; ih++) {
        const HitData_t* hit = mc->fListOfHits[ih];
        if (hit->fDeltaIndex >= 0) mc->fNHitsDelta += 1;
      }

      if (mc->fPdgID == PDGCode::e_minus) {
        float mom = mc->Momentum();
        if (mom < 20) {
          _nDeltaHitsTot  += mc->NHits();
          _nDeltaHitsReco += mc->fNHitsDelta;
        }
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
// Mode = 1:
// Mode = 2: event
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::debug(void* Data, int Mode) {
//-----------------------------------------------------------------------------
// print DeltaSeeds - pieces of delta electrons reconstructed within one station
//-----------------------------------------------------------------------------
    _data = (Data_t*) Data;

    if (Mode != 2) return -1; // beginRun not handled yet

    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
        _eventNumber       = en;
        int rc = InitMcDiag();
        if (rc < 0) return -1;
        associateMcTruth();
      }
    }

    if (_printComboHits) {
//-----------------------------------------------------------------------------
// print combo hits Z-ordered
//-----------------------------------------------------------------------------
      for (int is=0; is<kNStations; is++) {
        printHitData(nullptr,1);
        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = &_data->fFaceData[is][face];
          int nhits = fz->fHitData.size();
          for (int ih=0; ih<nhits; ih++) {
            HitData_t* hd = &fz->fHitData[ih];
            printHitData(hd,-1);
          }
        }
      }

      // HlPrint* hlp = HlPrint::Instance();
      // hlp->SetEvent(_data->event);
      // hlp->printComboHitCollection(_data->chCollTag.encode().data(),
      //                              _data->chfCollTag.encode().data(),
      //                              _data->sdmcCollTag.encode().data());
    }

    if (_printDeltaSeeds != 0) {
      for (int st=0; st<kNStations; ++st) {
        int nseeds = _data->NSeeds(st);
        printf("* station: %2i N(seeds): %3i\n",st,nseeds);
        if (nseeds > 0) {

          printDeltaSeed(nullptr,nullptr,1);

          for (int iseed=0; iseed<nseeds; ++iseed) {
            DeltaSeed* seed = _data->deltaSeed(st,iseed);
            DeltaSeedPar_t* sp = &_seedPar[st].at(iseed);
            printDeltaSeed(seed,sp,-1);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// print reconstructed delta candidates
//-----------------------------------------------------------------------------
    if (_printDeltaCandidates != 0) {
      int nd = _data->nDeltaCandidates();
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* [DeltaFinder::debug] N(delta candidates) = %5i\n",nd);
      printf("* ----------------------------------------------------------------------------------------------------------------\n");

      if (nd > 0) {
        for (int i=0; i<nd; i++) {
          DeltaCandidate* dc = _data->deltaCandidate(i);
          int pdg_id = -1;
          int sim_id = -1;
          if (dc->fMcPart) {
            pdg_id = dc->fMcPart->fPdgID;
            sim_id = dc->fMcPart->fID;
          }
          printf("* :dc:%05i  nh n(CE) ns s1  s2     X        Y        Z     chi21   chi22   htmin   htmax   t0    PdgID N(MC hits) simID=%5i pdgID=%10i \n",
                 dc->Index(),sim_id,pdg_id);
          printf("----------------------------------------------------------------------------------------------------------------\n");
          printf("            %3i  %3i %3i",dc->fNHits,dc->fNHitsCE,dc->fNSeeds);
          printf(" %2i  %2i %7.2f %7.2f %9.2f",dc->fFirstStation,dc->fLastStation,
                 dc->CofM.x(),dc->CofM.y(),dc->CofM.z());
          printf("                %30i %5i",pdg_id,dc->fNHitsMcP);
          printf("\n");
          printf("----------------------------------------------------------------------------------------------------------------\n");

          for (int is=dc->fFirstStation;is<=dc->fLastStation; is++) {
            DeltaSeed* ds = dc->seed[is];
            if (ds != NULL) {

              DeltaSeedPar_t* sp = &_seedPar[is].at(ds->Index());

              int face0 = ds->SFace(0);
              int face1 = ds->SFace(1);

              const HitData_t* hd0 = ds->HitData(face0);
              const HitData_t* hd1 = (face1 >= 0) ? ds->HitData(face1) : nullptr;

              printf("            %3i  %3i    %3i:%03i",ds->fNHits,sp->fNHitsCE,is,ds->Index());
              printf(" %7.2f %7.2f %9.2f",ds->CofM.x(),ds->CofM.y(),ds->CofM.z());
              float chi22 = (hd1) ? hd1->fChi2Min : -1;
              printf(" %7.1f %7.1f",hd0->fChi2Min, chi22);
              printf(" %7.1f %7.1f",ds->MinHitTime(),ds->MaxHitTime());
              printf(" %7.1f ",dc->T0(is));

              // printf("(");
              for (int face=0; face<kNFaces; face++) {
                const HitData_t* hd = ds->HitData(face);
                if (hd == nullptr) printf("(%5i:%5i)",-1,-1);
                else {
                  const ComboHit* hit = hd->fHit;
                  McPart_t* mcp = sp->fMcPart[face];
                  int mcid = mcp->fID;
                  printf("(%5i:%5i)",hit->strawId().asUint16(),mcid);
                }
              }
              printf("\n");
            }
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// print MC electrons
//-----------------------------------------------------------------------------
    if (_printElectrons) {
      int nmc = _list_of_mc_particles.GetEntriesFast();

      for (int i=0; i<nmc; i++) {
        McPart_t* mc = (McPart_t*) _list_of_mc_particles.At(i);

        if (((mc->fPdgID     == PDGCode::e_minus      ) or (mc->fPdgID     == PDGCode::e_plus)) &&
            (mc->Momentum() >  _printElectronsMinMom  ) &&
            (mc->Momentum() <  _printElectronsMaxMom  ) &&
            (mc->NHits()    >= _printElectronsMinNHits)    ) {

          float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

          if (fr < _printElectronsMaxFReco) {

            int delta_id(-1), nseeds(0);
            if (mc->fDelta) {
              delta_id = mc->fDelta->Index();
              nseeds   = mc->fDelta->NSeeds();
            }

            printf("* event: %4i electron.sim.id: %10i",_data->event->event(),mc->fID);
            printf(" mom = %7.3f time: %8.3f deltaID: %3i nseeds: %2i nhits: %3i/%3i stations:%2i:%2i",
                   mc->Momentum(), mc->Time(),
                   delta_id, nseeds,
                   mc->fNHitsDelta,
                   mc->NHits(),
                   mc->fFirstStation, mc->fLastStation);
            printf(" freco: %6.3f\n",fr);

            if (_printElectronsHits > 0) {
              int nh = mc->fListOfHits.size();
              if (nh > 0) printHitData(NULL,1);
              for (int ih=0; ih<nh; ih++) {
                printHitData(mc->fListOfHits[ih],-1);
              }
            }
          }
        }
      }
    }

    if (_printOTracker) printOTracker();

    if (_printShcol   ) printComboHitCollection();

    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printHitData(const HitData_t* Hd, int PrintBanner) {

    if (PrintBanner >= 0) {
      printf("#----------------------------------------------------------------------------------------------");
      printf("--------------------------------------------------------------------------------------\n");
      printf("#   I  SHID    flag   nsh St:F:Pl:Pn:Str   Time     TCorr     dt     eDep     wdist     wres   ");
      printf("     PDG simID    mom_PDG mom_ID     p      X        Y         Z   DeltaID   rad  edep\n");
      printf("#----------------------------------------------------------------------------------------------");
      printf("--------------------------------------------------------------------------------------\n");
      if (PrintBanner > 0) return;
    }

    const ComboHit* ch0 = &_data->chcol->at(0);
    const ComboHit* ch  = Hd->fHit;
    int loc             = ch-ch0;

    const StrawHitFlag* flag = &(*_data->outputChfColl)[loc];

    int radselOK        = (! flag->hasAnyProperty(StrawHitFlag::radsel   ));
    int edepOK          = (! flag->hasAnyProperty(StrawHitFlag::energysel));

    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999), mom_id(-9999), mom_pdg_id(-9999);
    float              mc_mom(-9999.);

    if (_mcDiag) {
      sim        = _mcUtils->getSimParticle(_data->event,ch->index(0));
      pdg_id     = _mcUtils->getPdgID(sim);
      sim_id     = _mcUtils->getID(sim);
      mc_mom     = _mcUtils->getStartMom(sim);
      mom_id     = _mcUtils->getMotherID(sim);
      mom_pdg_id = _mcUtils->getMotherPdgID(sim);
    }

    printf("%5i",loc);
    printf(" %5i 0x%08x %2i" ,ch->strawId().asUint16(),*((int*) flag),ch->nStrawHits());

    printf(" %2i:%i:%2i %2i %2i %8.2f %8.2f %7.3f %8.5f %8.3f %8.3f %10i %5i %10i %5i %8.3f %8.3f %8.3f %9.3f %5i %5i %5i\n",
           ch->strawId().station(),
           Hd->fZFace,
           ch->strawId().plane(),
           ch->strawId().panel(),
           ch->strawId().straw(),
           ch->time(),
           ch->correctedTime(),
           -1.,                             //sh->dt(),// ** FIXME!
           ch->energyDep(),
           ch->wireDist(),
           ch->posRes(ComboHit::wire),
           pdg_id,
           sim_id,
           mom_pdg_id,
           mom_id,
           mc_mom,
           ch->pos().x(),
           ch->pos().y(),
           ch->pos().z(),
           Hd->fDeltaIndex,
           radselOK,
           edepOK);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printOTracker() {

    int nhitso = 0;
    for (int is=0; is<kNStations; is++) {
      for (int face=0; face<kNFaces; face++) {
        FaceZ_t* fz = &_data->fFaceData[is][face];
        printf("#        --------------- station: %2i face: %2i nhits:%3li\n",
               is,face,fz->fHitData.size());
        int nh = fz->fHitData.size();
        if (nh > 0) printHitData(NULL,-1);
        for (int ih=0; ih<nh; ih++) {
          printHitData(&fz->fHitData[ih],face);
        }
        nhitso += nh;
      }
    }

    printf(" nhits, nhitso : %6i %6i \n", (int) _data->chcol->size(),nhitso);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printComboHit(const ComboHit* Ch, int PrintBanner) {

    if (PrintBanner >= 0) {
      printf("#-----------------------------------------------------------------------------");
      printf("----------------------------------------------------------\n");
      printf("#  I   SHID    flag   Plane Panel Straw    Time    T(corr)       dt      eDep ");
      printf("        PDG    ID       mPDG     mID       p  radsel  edep\n");
      printf("#-----------------------------------------------------------------------------");
      printf("----------------------------------------------------------\n");
      if (PrintBanner > 0) return;
    }

    const ComboHit* ch0 = &_data->chcol->at(0);
    int loc             = Ch-ch0;

    const StrawHitFlag* shf = &_data->outputChfColl->at(loc);

    int radsel = shf->hasAnyProperty(StrawHitFlag::radsel);
    int edep   = shf->hasAnyProperty(StrawHitFlag::energysel);

    const SimParticle* sim(nullptr);
    int                pdg_id(-9999), sim_id(-9999), mom_id(-9999), mom_pdg_id(-9999);
    float              mc_mom(-9999.);

    if (_mcDiag) {
      sim        = _mcUtils->getSimParticle(_data->event,Ch->index(0));
      pdg_id     = _mcUtils->getPdgID(sim);
      sim_id     = _mcUtils->getID(sim);
      mc_mom     = _mcUtils->getStartMom(sim);
      mom_id     = _mcUtils->getMotherID(sim);
      mom_pdg_id = _mcUtils->getMotherPdgID(sim);
    }

    printf("%5i",loc);
    printf(" %5i 0x%08x" ,Ch->strawId().asUint16(),*((int*) shf));

    printf(" %4i %4i %5i  %8.3f  %8.3f  %8.3f  %8.5f %10i %5i %10i    %5i %8.3f %5i %5i\n",
           Ch->strawId().plane(),
           Ch->strawId().panel(),
           Ch->strawId().straw(),
           Ch->time(),
           Ch->correctedTime(),
           -1.,                //Sh->dt(),//FIXME!
           Ch->energyDep(),
           pdg_id,
           sim_id,
           mom_pdg_id,
           mom_id,
           mc_mom,
           radsel,edep);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printComboHitCollection() {

    int nh = _data->chcol->size();
    printf("[DeltaFinderDiag::printComboHitCollection] : nhits : %6i\n", nh);

    printComboHit(nullptr,1);
    for (int i=0; i<nh; i++) {
      const ComboHit* ch = &_data->chcol->at(i);
      printComboHit(ch,-1);
    }

  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printDeltaSeed(DeltaSeed* seed, DeltaSeedPar_t* Sp, int PrintBanner) {

    if (PrintBanner >= 0) {
      printf("* -----------------------------------------------------------------------------------------");
      printf("-------------------------------------------------------------------------------\n");
      printf("* st seed  good type delta  SHID: MCID   SHID: MCID   SHID: MCID   SHID: MCID");
      printf("  chi2all/N chi2par/N chi2perp/N chi21   chi22 mintime  maxtime  <edep>      ");
      printf("X        Y         Z   nch nsh\n");
      printf("* -------------------------------------------------------------------------------------------");
      printf("-------------------------------------------------------------------------------\n");
    }

    if (PrintBanner > 0) return;

    printf("* %2i  %03i %5i   %02i",seed->Station(),seed->Index(),seed->fGood,seed->fType);
    printf(" %5i",seed->fDeltaIndex);
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
    for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd = seed->HitData(face);
      if (hd == nullptr) printf("(%5i:%5i)",-1,-1);
      else {
        const ComboHit* hit = hd->fHit;
        McPart_t* mcp = Sp->fMcPart[face];
        int mcid      = mcp->fID;
        printf("(%5i:%5i)",hit->strawId().asUint16(),mcid);
      }
    }

    printf(" %8.2f   %7.2f  %7.2f  %7.2f %7.2f",seed->Chi2TotN(),seed->Chi2ParN(),seed->Chi2PerpN(),seed->fChi21,seed->fChi22);
    printf("%8.1f %8.1f %8.5f",seed->MinHitTime(),seed->MaxHitTime(),seed->EDep());
    printf(" %8.3f %8.3f %9.3f",seed->CofM.x(),seed->CofM.y(),seed->CofM.z());
    //    printf("%4i",seed->fNFacesWithHits);
    printf("%4i",seed->NHits());
    printf("%4i",seed->NStrawHits());
    printf("\n");
//-----------------------------------------------------------------------------
// after that, print combo hits
//-----------------------------------------------------------------------------
    printHitData(nullptr,1);
    for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd = seed->HitData(face);
      if (hd == nullptr)                                              continue;
      printHitData(hd,-1);
    }

  }

}

DEFINE_ART_CLASS_TOOL(mu2e::DeltaFinderDiag)

#endif
