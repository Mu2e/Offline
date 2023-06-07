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
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CalPatRec/inc/HlPrint.hh"
#include "Offline/CalPatRec/inc/DeltaFinderAlg.hh"
#include "Offline/CalPatRec/inc/McPart_t.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

using namespace std;

namespace mu2e {

  using namespace DeltaFinderTypes;

  class SimParticle;
  class StrawDigiMCCollection;

  class DeltaFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets  =  10,
      kNSeedHistSets   = 100,
      kNSeed2HistSets  = 100,
      kNDeltaHistSets  = 100,
      kNProtonHistSets = 100,
      kNMcHistSets     = 500
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
      TH1F*  fChi2N;
      TH1F*  fChi2ParN;
      TH1F*  fChi2PerpN;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNSeeds;
      TH2F*  fNSeedsVsStation;
      TH1F*  fNMc;
      TH1F*  fPDGCode;
      TH1F*  fNMcHits;
      TH1F*  fNDelta;
      TH1F*  fNDeltaHitsT;
      TH1F*  fNDeltaHitsR;
      TH2F*  fNPRecoVsMc;
      TH2F*  fNP10RecoVsMc;
      TH2F*  fNP15RecoVsMc;
      TH2F*  fNP10RecoVsPpi;
      TH2F*  fNP15RecoVsPpi;
      TH1F*  fNch;
      TH1F*  fNChFlaggedDelta;
      TH1F*  fNChFlaggedProton;
    };

    struct SeedHist_t {
      TH1F*  fChi2TotN;
      TH1F*  fChi2ParN;
      TH1F*  fChi2PerpN;
      TH1F*  fChi2Delta;
      TH1F*  fChi2DeltaPar;
      TH1F*  fChi2DeltaPerp;
      TH1F*  fChi2SeedHit;                        // chi2 of the first two hits (along the wire)
      TH1F*  fChi2Neighbour;
      TH1F*  fChi2Radial;
      TH1F*  fNHitsPerFace;
      TH1F*  fNHitsPerSeed;
      TH1F*  fSeedRadius;
      TH1F*  fSeedMomentum;
      TH1F*  fNHits;
      TH1F*  fNHitsCE;
      TH1F*  fNSim;
      TH1F*  fNMom;
      TH1F*  fEDep;
      TH1F*  fDtDelta;
      TH1F*  fDt12;
      TH1F*  fDtCorr12;

      TH1F*  fChi2H;
      TH1F*  fChi2HPar;
      TH1F*  fChi2HPerp;
    };

    struct Seed2Hist_t {
      TH1F*  fDt;                       // dt between the two seeds in teh same station
    };

    struct ProtonHist_t {
      TH1F*  fNComboHits;
      TH1F*  fNStrawHits;
      TH1F*  fNHitsCE;
      TH1F*  fNHitsF;
      TH1F*  fNStationsWithHits;
      TH1F*  fMcMom;
      TH1F*  fPDGCode;
      TH1F*  fFFlagged;
      TH1F*  fDPhi;
    };
                                        // hits referred to here are the combo hits
    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom;
      TH1F*  fNHits;
      TH1F*  fNHitsDelta;
      TH1F*  fNHitsProton;
      TH1F*  fNChFlaggedDelta;
      TH1F*  fNChFlaggedProton;
      TH2F*  fNHfdVsNH;
      TH2F*  fNHfpVsNH;
      TH1F*  fFractReco;
      TH1F*  fMaxSeg;
      TH2F*  fFractRecoVsNHits;
      TH1F*  fHitDt;
    };

    struct Hist_t {
      EventHist_t*  fEvent [kNEventHistSets ];
      DeltaHist_t*  fDelta [kNDeltaHistSets ];
      McHist_t*     fMc    [kNMcHistSets    ];
      ProtonHist_t* fProton[kNProtonHistSets];
      SeedHist_t*   fSeed  [kNSeedHistSets  ];
      Seed2Hist_t*  fSeed2 [kNSeed2HistSets ];
    };
//-----------------------------------------------------------------------------
// additional parameters of the DeltaCandidate , for analysis
//-----------------------------------------------------------------------------
    struct DeltaCandidatePar_t {
      McPart_t*    fMcPart;               // "best" MC particle
      int          fNHitsMcP;             // N combo hits by the "best" particle"
      int          fNHitsCE;              // N(hits) by CE
      float        dxy[kNStations];       // coordinate residuals ... TBD
      float        fChi2N;                // all-hit delta chi2s
      float        fChi2ParN;             // parallel component
      float        fChi2PerpN;            // perpendicular component
    };

    std::vector<DeltaCandidatePar_t>  _deltaPar;
//-----------------------------------------------------------------------------
// additional parameters of the ProtonCandidate, for analysis
//-----------------------------------------------------------------------------
    struct ProtonCandidatePar_t {
      McPart_t*    fMcPart;               // "best" MC particle
      int          fNHitsMcP;             // N combo hits by the "best" particle"
      int          fNHitsCE;              // N(hits) produced by CE
      int          fNFlaggedHits;         // N(hits) flagged as proton

      std::vector<McPart_t*> fHitMcPart[kNStations][kNFaces]; // reset properly !
    };

    std::vector<ProtonCandidatePar_t>  _protonPar;
    TObjArray                          _listOfMcProtons;
//-----------------------------------------------------------------------------
// additional parameters of the DeltaSeed , for analysis
//-----------------------------------------------------------------------------
    struct DeltaSeedPar_t {
      float        fDtDelta;              // time residual to delta
      float        fDt12;                 // time residual of the two hits
      float        fDtCorr12;             // time residual of the two hits
      McPart_t*    fMcPart[kNFaces];
      int          fNSim;                 // N different simID's cof the hits
      int          fNMom;                 // N(mothers) of the particles which produced the hits
      int          fNHitsCE;              // number of associated CE hits
      McPart_t*    fPreSeedMcPart[2];     // McPart_t's for initial stereo intersection
      float        fDeltaChi2;
      float        fDeltaChi2Par;
      float        fDeltaChi2Perp;
                                          // hit chi2s wrt the segment
      float        fChi2H    [kNFaces];
      float        fChi2HPar [kNFaces];
      float        fChi2HPerp[kNFaces];

      int MCTruth () { return (fPreSeedMcPart[0] != NULL) && (fPreSeedMcPart[0] == fPreSeedMcPart[1]) ; }
    };

    std::vector<DeltaSeedPar_t>           _seedPar[kNStations];

  protected:

    bool           _mcDiag;
    int            _printOTracker;
    int            _printComboHits;
    int            _printGoodComboHits;
    int            _printElectrons;
    int            _printElectronsHits;
    int            _printElectronsMinNHits;
    float          _printElectronsMaxFReco;
    float          _printElectronsMinMom;
    float          _printElectronsMaxMom;
    int            _printDeltaSeeds;
    int            _printDeltaCandidates;
    int            _printShcol;
    int            _printSeedNParents;   // n(hits) for seeds with > 2 parents

    int            _printMcProtons;
    int            _printProtonHits;
    int            _printProtonSeeds;
    int            _printProtonCandidates;

    int            _eventNumber;
    int            _nDeltaHitsTot;
    int            _nDeltaHitsReco;

    TObjArray      _listOfMcParticles; // list of particles with at least one ComboHit in the tracker
    TObjArray      _listOfMcPartHit ;  // for each ComboHit, pointer to its McPart

    int            _nMcProtons;
    int            _nMcProtons10;
    int            _nMcProtons15;
    int            _nRecoProtons;
    int            _nRecoProtons10;
    int            _nRecoProtons15;

    int            _nChFlaggedDelta;
    int            _nChFlaggedProton;

    Data_t*        _data;                 // shared data, passed from the caller
    Hist_t         _hist;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    float          _ppi;                // proton pulse intensity

  public:

    DeltaFinderDiag(const fhicl::Table<mu2e::DeltaFinderTypes::Config>& config);
    ~DeltaFinderDiag();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
    McPart_t*   mcPart   (int I) { return (McPart_t*) _listOfMcParticles.At(I); }
    McPart_t*   mcPartHit(int I) { return (McPart_t*) _listOfMcPartHit.At  (I); }
    McPart_t*   mcProton (int I) { return (McPart_t*) _listOfMcProtons.At  (I); }
//-----------------------------------------------------------------------------
// other functions
//-----------------------------------------------------------------------------
    void        bookDeltaHistograms (DeltaHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookEventHistograms (EventHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookMcHistograms    (McHist_t*     Hist, art::TFileDirectory* Dir);
    void        bookSeedHistograms  (SeedHist_t*   Hist, art::TFileDirectory* Dir);
    void        bookSeed2Histograms (Seed2Hist_t*  Hist, art::TFileDirectory* Dir);
    void        bookProtonHistograms(ProtonHist_t* Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms (EventHist_t*  Hist);
    void        fillDeltaHistograms (DeltaHist_t*  Hist, DeltaCandidate*  Delta, DeltaCandidatePar_t*  Dcp);
    void        fillMcHistograms    (McHist_t*     Hist, McPart_t*        Mc   );
    void        fillSeedHistograms  (SeedHist_t*   Hist, DeltaSeed*       Seed , DeltaSeedPar_t* SeedPar);
    void        fillSeed2Histograms (Seed2Hist_t*  Hist, DeltaSeed*       Seed1, DeltaSeed* Seed2  );
    void        fillProtonHistograms(ProtonHist_t* Hist, ProtonCandidate* Prot , ProtonCandidatePar_t* Pcp);

    McPart_t*   findParticle (const SimParticle* Sim);

    int         associateMcTruth();
    int         initDeltaCandMcParameters ();
    int         initDeltaSeedMcParameters ();
    int         initProtonCandMcParameters();
    int         initMcDiag                ();

    int         precalculateRecoParameters();

    void        printDeltaCandidates();

    void        printSegment (DeltaSeed* Seed, DeltaSeedPar_t* Sp, int PrintBanner);
    void        printHitData (const HitData_t* Hd, int PrintBanner);
    void        printOTracker();

    void        printComboHit (const ComboHit* Sh, int PrintBanner);
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
  DeltaFinderDiag::DeltaFinderDiag(const fhicl::Table<mu2e::DeltaFinderTypes::Config>& config):
    _mcDiag                (config().mcDiag()                ),
    _printOTracker         (config().printOTracker()         ),
    _printComboHits        (config().printComboHits()        ),
    _printGoodComboHits    (config().printGoodComboHits()    ),
    _printElectrons        (config().printElectrons()        ),
    _printElectronsHits    (config().printElectronsHits()    ),
    _printElectronsMinNHits(config().printElectronsMinNHits()),
    _printElectronsMaxFReco(config().printElectronsMaxFReco()),
    _printElectronsMinMom  (config().printElectronsMinMom()  ),
    _printElectronsMaxMom  (config().printElectronsMaxMom()  ),
    _printDeltaSeeds       (config().printDeltaSeeds()       ),
    _printDeltaCandidates  (config().printDeltaCandidates()  ),
    _printShcol            (config().printShcol()            ),
    _printSeedNParents     (config().printSeedNParents()     ),
    _printMcProtons        (config().printMcProtons()        ),
    _printProtonHits       (config().printProtonHits()       ),
    _printProtonSeeds      (config().printProtonSeeds()      ),
    _printProtonCandidates (config().printProtonCandidates() )
  {
    printf(" DeltaFinderDiag::DeltaFinderDiag : HOORAY Config! \n");

    if (_mcDiag != 0) _mcUtils = art::make_tool  <McUtilsToolBase>(config().mcUtils,"mcUtils");
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  DeltaFinderDiag::~DeltaFinderDiag() {
  }

//-----------------------------------------------------------------------------
// for each DeltaSeed, create a list of SimParticle*'s parallel to its list of combo hits
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::associateMcTruth() {

    if (_data->chcol->size() == 0) return 0;

    initDeltaSeedMcParameters();

    initDeltaCandMcParameters();

    initProtonCandMcParameters();
//-----------------------------------------------------------------------------
// now, for each MC electron calculate the number of hits flagged as delta
//-----------------------------------------------------------------------------
    int nmc    = _listOfMcParticles.GetEntriesFast();

    _nDeltaHitsTot    = 0;
    _nDeltaHitsReco   = 0;

    _nChFlaggedProton = 0;
    _nChFlaggedDelta  = 0;

    const ComboHit* ch_hit_0 = &_data->chcol->at(0);

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _listOfMcParticles.At(i);
      mc->fNHitsDelta       = 0;
      mc->fNHitsProton      = 0;
      mc->fNChFlaggedDelta  = 0;
      mc->fNChFlaggedProton = 0;
//-----------------------------------------------------------------------------
// loop over the hits of MC delta electron and calculate fraction of them which have
// been reconstructed as hits of reconstructed delta electrons
//-----------------------------------------------------------------------------
      int nh = mc->fListOfHits.size();
      for (int ih=0; ih<nh; ih++) {
        const HitData_t* hd = mc->fListOfHits[ih];
        int loc = hd->fHit-ch_hit_0;
        const StrawHitFlag* flag = &_data->outputChColl->at(loc).flag();

        if (flag->hasAnyProperty(StrawHitFlag::bkg)) mc->fNChFlaggedDelta += 1;

        assert ((hd->fDeltaIndex >= 0) == flag->hasAnyProperty(StrawHitFlag::bkg));

        if (! flag->hasAnyProperty(StrawHitFlag::energysel)) mc->fNChFlaggedProton += 1;

        if (hd->fProtonIndex > 0) {
          ProtonCandidate* pc = _data->protonCandidate(hd->fProtonIndex);
//-----------------------------------------------------------------------------
// this is just a consistency check - 1-station proton candidates are required
// to have eDep > 4 keV
//-----------------------------------------------------------------------------
          if ((pc->nStationsWithHits() == 1) and(pc->eDep() > 0.004)) {
            assert (not flag->hasAnyProperty(StrawHitFlag::energysel));
          }
        }
      }

      _nChFlaggedDelta  += mc->fNChFlaggedDelta;
      _nChFlaggedProton += mc->fNChFlaggedProton;

      if (mc->fDelta ) mc->fNHitsDelta  = mc->fDelta->nHits();
      if (mc->fProton) mc->fNHitsProton = mc->fProton->nHitsTot();

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
  McPart_t* DeltaFinderDiag::findParticle(const SimParticle* Sim) {
    McPart_t* found(0);

    int n = _listOfMcParticles.GetEntriesFast();

    for (int i=0; i<n; i++) {
      McPart_t* mc = (McPart_t*) _listOfMcParticles.At(i);
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
    Hist->fChi2N      = Dir->make<TH1F>("chi2n"     , "Chi2/N"          , 1000, 0., 100.);
    Hist->fChi2ParN   = Dir->make<TH1F>("chi2parn"  , "Chi2Par/N"       , 1000, 0., 100.);
    Hist->fChi2PerpN  = Dir->make<TH1F>("chi2perpn" , "Chi2Perp/N"      , 1000, 0., 100.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number", 100, 0., 100000.);
    Hist->fNSeeds          = Dir->make<TH1F>("nseeds", "N(seeds)"   , 200, 0., 2000.);
    Hist->fNSeedsVsStation = Dir->make<TH2F>("ns_vs_st", "N(seeds) vs station", 20, 0., 20.,200,0,200);

    Hist->fNMc             = Dir->make<TH1F>("nmc"       , "N(MC particles)", 100, 0., 1000.);
    Hist->fPDGCode         = Dir->make<TH1F>("pdg_code"  , "PDG Code"       ,1000, -2500., 2500.);
    Hist->fNMcHits         = Dir->make<TH1F>("n_mc_hits" , "N(MC hits)"     , 250, 0., 500.);
    Hist->fNDelta          = Dir->make<TH1F>("n_delta"   , "N(reco deltas)" , 250, 0., 250.);
    Hist->fNDeltaHitsT     = Dir->make<TH1F>("n_delta_ht", "N(delta hits T)", 500, 0., 5000.);
    Hist->fNDeltaHitsR     = Dir->make<TH1F>("n_delta_hr", "N(delta hits R)", 500, 0., 5000.);

    Hist->fNPRecoVsMc    = Dir->make<TH2F>("np_reco_vs_mc"  , "Nprot: reco vs MC)"  , 100, 0., 100.,100,0,100);
    Hist->fNP10RecoVsMc  = Dir->make<TH2F>("np10_reco_vs_mc", "Np10: reco vs MC)", 100, 0., 100.,100,0,100);
    Hist->fNP15RecoVsMc  = Dir->make<TH2F>("np15_reco_vs_mc", "Np15: reco vs MC)", 100, 0., 100.,100,0,100);
    Hist->fNP10RecoVsPpi = Dir->make<TH2F>("np10_reco_vs_ppi", "Np10r vs PPI)", 20, 0., 20.,50,0,50);
    Hist->fNP15RecoVsPpi = Dir->make<TH2F>("np15_reco_vs_ppi", "Np15r vs PPI)", 20, 0., 20.,50,0,50);

    Hist->fNch              = Dir->make<TH1F>("nch"   , "N(CH) total"             , 500, 0., 10000.);
    Hist->fNChFlaggedDelta  = Dir->make<TH1F>("nch_fd", "N(CH) flagged as delta"  , 500, 0., 10000.);
    Hist->fNChFlaggedProton = Dir->make<TH1F>("nch_fp", "N(CH) flagged as proton)", 500, 0., 10000.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookMcHistograms(McHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fPDGCode     = Dir->make<TH1F>("pdg"  , "PDG code"        ,1000, -2500, 2500);
    Hist->fMom         = Dir->make<TH1F>("mom"  , "momentum"        , 500, 0.,1000.);
    Hist->fNHits       = Dir->make<TH1F>("nhits", "N(hits)"         , 100, 0., 100.);
    Hist->fNHitsDelta  = Dir->make<TH1F>("nh_rd", "N(hits on a matching reconstructed delta)" , 100, 0., 100.);
    Hist->fNHitsProton = Dir->make<TH1F>("nh_rp", "N(hits on a matching reconstructed proton)", 100, 0., 100.);

    Hist->fNChFlaggedDelta  = Dir->make<TH1F>("nch_fd", "N(CH) flagged as delta"  , 100, 0., 100.);
    Hist->fNChFlaggedProton = Dir->make<TH1F>("nch_fp", "N(CH) flagged as proton)", 100, 0., 100.);

    Hist->fFractReco  = Dir->make<TH1F>("fractr", "NR/N"           , 100, 0.,   1.);
    Hist->fMaxSeg     = Dir->make<TH1F>("max_seg", "Max N Segments",  20, 0.,  20.);
    Hist->fHitDt      = Dir->make<TH1F>("hit_dt" , "Hit TMax-TMin" , 100, 0., 200.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);

    Hist->fNHfpVsNH = Dir->make<TH2F>("nhfp_vs_nh", "Nh flagged P vs nh", 100, 0., 100.,100,0,100);
    Hist->fNHfdVsNH = Dir->make<TH2F>("nhfd_vs_nh", "Nh flagged D vs nh", 100, 0., 100.,100,0,100);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookProtonHistograms(ProtonHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fNComboHits = Dir->make<TH1F>("nch"       , "N(combo hits)"   , 200, 0., 200.);
    Hist->fNStrawHits = Dir->make<TH1F>("nsh"       , "N(straw hits)"   , 200, 0., 200.);
    Hist->fNHitsCE    = Dir->make<TH1F>("nch_ce"    , "N(combo hits CE)",  50, 0.,  50.);
    Hist->fNHitsF     = Dir->make<TH1F>("nchf"      , "N(flagged CH)"   ,  50, 0.,  50.);
    Hist->fPDGCode    = Dir->make<TH1F>("pdg"       , "PDG code"        ,1000, -2500., 2500.);
    Hist->fMcMom      = Dir->make<TH1F>("mc_mom"    , "N(hits)"         , 500, 0., 1000.);
    Hist->fFFlagged   = Dir->make<TH1F>("fflagged"  , "N(flagged)/N"    , 100, 0.,    1);
    Hist->fDPhi       = Dir->make<TH1F>("dphi"      , "Delta Phi"       , 200, -1.,   1);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookSeedHistograms(SeedHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fChi2TotN       = Dir->make<TH1F>("chi2totn"       , "Chi2Tot/N"           , 1000, 0., 100.);
    Hist->fChi2ParN       = Dir->make<TH1F>("chi2parn"       , "Chi2Par/N"           , 1000, 0., 100.);
    Hist->fChi2PerpN      = Dir->make<TH1F>("chi2perpn"      , "Chi2Perp/N"          , 1000, 0., 100.);
    Hist->fChi2Delta      = Dir->make<TH1F>("chi2_delta"     , "Chi2Delta"           , 1000, 0., 100.);
    Hist->fChi2DeltaPar   = Dir->make<TH1F>("chi2_delta_par" , "Chi2DeltaPar"        , 1000, 0., 100.);
    Hist->fChi2DeltaPerp  = Dir->make<TH1F>("chi2_delta_perp", "Chi2DeltaPerp"       , 1000, 0., 100.);
    Hist->fChi2SeedHit    = Dir->make<TH1F>("chi2_seed_hit"  , "Chi seed stereo hit" , 1000, 0.,  50.);
    Hist->fChi2Neighbour  = Dir->make<TH1F>("chi2_nb"        , "Chi2 neighbour"      , 1000, 0.,  50.);
    Hist->fChi2Radial     = Dir->make<TH1F>("chi2_r"         , "Chi2 radial"         , 1000, 0.,  50.);
    Hist->fNHitsPerFace   = Dir->make<TH1F>("nhits_face"     , "Nhits per face"      , 20, 0., 20.);
    Hist->fNHits          = Dir->make<TH1F>("nhits"          , "Number of hits"      , 10, 0., 10.);
    Hist->fNHitsCE        = Dir->make<TH1F>("nhits_ce"       , "Number of CE hits"   , 10, 0., 10.);
    Hist->fSeedRadius     = Dir->make<TH1F>("rad"            , "Seed radius"         , 500, 0., 1000.);
    Hist->fSeedMomentum   = Dir->make<TH1F>("mom"            , "Seed momentum"       , 400, 0., 400.);
    Hist->fNSim           = Dir->make<TH1F>("nsim"           , "N(sim particls)"     , 10,    0.,  10.);
    Hist->fNMom           = Dir->make<TH1F>("nmom"           , "N(MC moms)"          , 10,    0.,  10.);
    Hist->fEDep           = Dir->make<TH1F>("edep"           , "mean E(Dep)"         ,100,    0.,   0.01);
    Hist->fDtDelta        = Dir->make<TH1F>("dt_delta"       , "T(seed)-T(delta)"    ,200, -100., 100.);
    Hist->fDt12           = Dir->make<TH1F>("dt12"           , "T(h0)-T(h1)"         ,200, -100., 100.);
    Hist->fDtCorr12       = Dir->make<TH1F>("dtcorr12"       , "Tcorr(h0)-Tcorr(h1)" ,200, -100., 100.);

    Hist->fChi2H          = Dir->make<TH1F>("chi2h"          , "Chi2 hit"             , 1000, 0., 100.);
    Hist->fChi2HPar       = Dir->make<TH1F>("chi2h_par"      , "Chi2 hit par"         , 1000, 0., 100.);
    Hist->fChi2HPerp      = Dir->make<TH1F>("chi2h_perp"     , "Chi2 hit Perp"        , 1000, 0., 100.);
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
    book_mc_histset[  3] = 1;                // electrons with last>first
    book_mc_histset[  4] = 1;                // electrons with last>first and 4+ hits
    book_mc_histset[  5] = 1;                // electrons with last>first, 5+ hits, and reco delta
    book_mc_histset[  6] = 1;                // electrons with last>first, 5+ hits, and p < 20

    book_mc_histset[100] = 1;                // electrons and 20 < p < 80 MeV/c
    book_mc_histset[101] = 1;                // electrons and p > 80 MeV/c

    book_mc_histset[200] = 1;                // all protons
    book_mc_histset[201] = 1;                // protons with N(hits) > 10
    book_mc_histset[202] = 1;                // protons with N(hits) > 15

    for (int i=0; i<kNMcHistSets; i++) {
      if (book_mc_histset[i] != 0) {
        sprintf(folder_name,"mc_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fMc[i] = new McHist_t;
        bookMcHistograms(_hist.fMc[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book proton histograms
//-----------------------------------------------------------------------------
    int book_proton_histset[kNProtonHistSets];
    for (int i=0; i<kNProtonHistSets; i++) book_proton_histset[i] = 0;

    book_proton_histset[ 0] = 1;                // all proton candidates
    book_proton_histset[ 1] = 1;                // proton candidates N(combo hits) >= 10
    book_proton_histset[ 2] = 1;                // proton candidates N(combo hits) >= 15

    for (int i=0; i<kNProtonHistSets; i++) {
      if (book_proton_histset[i] != 0) {
        sprintf(folder_name,"proton_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fProton[i] = new ProtonHist_t;
        bookProtonHistograms(_hist.fProton[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillDeltaHistograms(DeltaHist_t* Hist, DeltaCandidate* Delta, DeltaCandidatePar_t* Dcp) {
    int n_seeds = Delta->fNSeeds;
    Hist->fNSeeds->Fill(n_seeds);
    Hist->fNHits->Fill(Delta->fNHits);
    Hist->fNHitsCE->Fill(Dcp->fNHitsCE);
    Hist->fNStrawHits->Fill(Delta->fNStrawHits);

    float mom(-1), pdg_code(-1.e6);
    if (Dcp->fMcPart) {
      mom      = Dcp->fMcPart->Momentum();
      pdg_code = Dcp->fMcPart->fPdgID;
    }

    Hist->fPDGCode->Fill(pdg_code);
    Hist->fMcMom->Fill(mom);
    Hist->fEDep->Fill(Delta->EDep());
    Hist->fNBestVsN->Fill(Delta->nHits(),Dcp->fNHitsMcP);

    float fbest = float(Dcp->fNHitsMcP)/Delta->nHits();
    Hist->fFBest->Fill(fbest);

    for(int is=0; is<kNStations; ++is) {
      Hist->fDxy->Fill(Dcp->dxy[is]);
    }

    Hist->fChi2N ->Fill   (Dcp->fChi2N    );
    Hist->fChi2ParN ->Fill(Dcp->fChi2ParN );
    Hist->fChi2PerpN->Fill(Dcp->fChi2PerpN);
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();

    Hist->fEventNumber->Fill(event_number);
    Hist->fNSeeds->Fill(_data->nSeedsTot());

    for (int is=0; is<kNStations; ++is) {

      int nseeds = _data->NSeeds(is);
      Hist->fNSeedsVsStation->Fill(is,nseeds);
    }

    int ndelta = _data->nDeltaCandidates();
    Hist->fNDelta->Fill(ndelta);
//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _listOfMcParticles.GetEntriesFast();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _listOfMcParticles.At(i);
      int n_mc_hits = mc->fListOfHits.size();

      Hist->fNMcHits->Fill(n_mc_hits);
      Hist->fPDGCode->Fill(mc->fPdgID);
    }

    Hist->fNMc->Fill(nmc);
    Hist->fNDeltaHitsT->Fill(_nDeltaHitsTot);
    Hist->fNDeltaHitsR->Fill(_nDeltaHitsReco);

    Hist->fNPRecoVsMc  ->Fill(_nMcProtons  ,_nRecoProtons  );
    Hist->fNP10RecoVsMc->Fill(_nMcProtons10,_nRecoProtons10);
    Hist->fNP15RecoVsMc->Fill(_nMcProtons15,_nRecoProtons15);

    Hist->fNP10RecoVsPpi->Fill(_ppi/1.e7,_nRecoProtons10);
    Hist->fNP15RecoVsPpi->Fill(_ppi/1.e7,_nRecoProtons15);

    int nch = _data->chcol->size();
    Hist->fNch->Fill(nch);
    Hist->fNChFlaggedDelta->Fill(_nChFlaggedDelta);
    Hist->fNChFlaggedProton->Fill(_nChFlaggedProton);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillMcHistograms(McHist_t* Hist, McPart_t* Mc) {

    float mom = Mc->Momentum();

    Hist->fPDGCode->Fill(Mc->fPdgID);
    Hist->fMom->Fill(mom);
    Hist->fNHits->Fill(Mc->NHits());

    Hist->fNHitsDelta->Fill(Mc->fNHitsDelta);
    Hist->fNHitsProton->Fill(Mc->fNHitsProton);

    Hist->fNChFlaggedDelta->Fill(Mc->fNChFlaggedDelta);
    Hist->fNChFlaggedProton->Fill(Mc->fNChFlaggedProton);

    Hist->fNHfdVsNH->Fill(Mc->NHits(),Mc->fNChFlaggedDelta);
    Hist->fNHfpVsNH->Fill(Mc->NHits(),Mc->fNChFlaggedProton);

    float freco = Mc->fNHitsDelta/(Mc->NHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNHits->Fill(Mc->NHits(),freco);

    int max_nseg = Mc->fLastStation-Mc->fFirstStation+1;
    Hist->fMaxSeg->Fill(max_nseg);

    float dt = Mc->HitDt();
    Hist->fHitDt->Fill(dt);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillProtonHistograms(ProtonHist_t* Hist, ProtonCandidate* Pc, ProtonCandidatePar_t* Pcp) {

    Hist->fNComboHits->Fill(Pc->nHitsTot());
    Hist->fNStrawHits->Fill(Pc->nStrawHitsTot());
    Hist->fNHitsCE->Fill(Pcp->fNHitsCE);
    Hist->fNHitsF->Fill(Pcp->fNFlaggedHits);

    float mom(-1), pdg_code(-1.e6);
    if (Pcp->fMcPart) {
      mom      = Pcp->fMcPart->Momentum();
      pdg_code = Pcp->fMcPart->fPdgID;
    }

    Hist->fPDGCode->Fill(pdg_code);
    Hist->fMcMom->Fill(mom);

    float f_flagged = float(Pcp->fNFlaggedHits)/Pc->nHitsTot();
    Hist->fFFlagged->Fill(f_flagged);

    for (int is=Pc->fFirstStation+1; is<=Pc->fLastStation; is++) {
      if ((Pc->nHitsStation(is-1) > 0) and (Pc->nHitsStation(is) > 0)) {
        float dphi = Pc->phi(is)-Pc->phi(is-1);
        Hist->fDPhi->Fill(dphi);
      }
    }
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillSeedHistograms(SeedHist_t* Hist, DeltaSeed* Seed, DeltaSeedPar_t* Par) {

    Intersection_t res;

    Hist->fChi2TotN ->Fill(Seed->Chi2TotN ());
    Hist->fChi2PerpN->Fill(Seed->Chi2PerpN());
    Hist->fChi2ParN ->Fill(Seed->Chi2ParN ());

    Hist->fChi2Delta    ->Fill(Par->fDeltaChi2    );
    Hist->fChi2DeltaPar ->Fill(Par->fDeltaChi2Par );
    Hist->fChi2DeltaPerp->Fill(Par->fDeltaChi2Perp);

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

      Hist->fChi2SeedHit->Fill(chi2_1);
      Hist->fChi2SeedHit->Fill(chi2_2);

      for (int face=0; face<kNFaces; face++) {
        const HitData_t* hd  = Seed->HitData(face);
        if (hd != nullptr) {
          Hist->fChi2H    ->Fill(Par->fChi2H    [face]);
          Hist->fChi2HPar ->Fill(Par->fChi2HPar [face]);
          Hist->fChi2HPerp->Fill(Par->fChi2HPerp[face]);
        }
      }
    }

    Hist->fSeedRadius->Fill    (Seed->CofM.rho());

    double mom (-1.);
    if (Par->fPreSeedMcPart[0]) mom = Par->fPreSeedMcPart[0]->Momentum();

    Hist->fSeedMomentum->Fill(mom);

    Hist->fNHits->Fill(Seed->nHits());
    Hist->fNHitsCE->Fill(Par->fNHitsCE);
    Hist->fNSim->Fill(Par->fNSim);
    Hist->fNMom->Fill(Par->fNMom);

    Hist->fEDep->Fill(Seed->EDep());

    Hist->fDtDelta->Fill(Par->fDtDelta);
    Hist->fDt12->Fill(Par->fDt12);
    Hist->fDtCorr12->Fill(Par->fDtCorr12);
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillSeed2Histograms(Seed2Hist_t* Hist, DeltaSeed* Seed1, DeltaSeed* Seed2) {
    float dt = Seed1->TMean()-Seed2->TMean();
    Hist->fDt->Fill(dt);
  }

//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;
//-----------------------------------------------------------------------------
// precalculate reco parameters - do it *before* precalculating MC parameters
//-----------------------------------------------------------------------------
    precalculateRecoParameters();
//-----------------------------------------------------------------------------
// precalculate MC-specific info
//-----------------------------------------------------------------------------
    art::Handle<mu2e::ProtonBunchIntensity> ppiH;
    _data->event->getByLabel("PBISim", ppiH);

    if (!ppiH.isValid()) {
      mf::LogWarning("DeltaFinderDiag::fillHistograms") << "no ProtonBunchIntensity objects found";
      _ppi = -1;
    }
    else {
      _ppi = ppiH->intensity();
    }

    int en = _data->event->event();
    if (_mcDiag) {
      if (_eventNumber != en) {
        _eventNumber = en;
        initMcDiag();
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

        fillSeedHistograms(_hist.fSeed[0],seed,sp);

        if      (seed->nHits() == 1) fillSeedHistograms(_hist.fSeed[ 1],seed,sp);
        else if (seed->nHits() == 2) fillSeedHistograms(_hist.fSeed[ 2],seed,sp);
        else if (seed->nHits() == 3) fillSeedHistograms(_hist.fSeed[ 3],seed,sp);
        else if (seed->nHits() == 4) fillSeedHistograms(_hist.fSeed[ 4],seed,sp);

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
            printSegment(seed,sp,0);
          }
        }
        DeltaCandidate* dc(nullptr);
        if ((seed->deltaIndex() >= 0) and (seed->fDeltaIndex < 10000)) {
          dc = _data->deltaCandidate(seed->deltaIndex());
        }
        if (dc) {
//-----------------------------------------------------------------------------
// this seed has been associated with a delta
//-----------------------------------------------------------------------------
          if (dc->nHits() >= 5) {
            if      (seed->nHits() == 1) fillSeedHistograms(_hist.fSeed[21],seed,sp);
            else if (seed->nHits() == 2) fillSeedHistograms(_hist.fSeed[22],seed,sp);
            else if (seed->nHits() == 3) fillSeedHistograms(_hist.fSeed[23],seed,sp);
            else if (seed->nHits() == 4) fillSeedHistograms(_hist.fSeed[24],seed,sp);
          }

          if (dc->NSeeds() > 3) {
            if      (seed->nHits() == 1) fillSeedHistograms(_hist.fSeed[31],seed,sp);
            else if (seed->nHits() == 2) fillSeedHistograms(_hist.fSeed[32],seed,sp);
            else if (seed->nHits() == 3) fillSeedHistograms(_hist.fSeed[33],seed,sp);
            else if (seed->nHits() == 4) fillSeedHistograms(_hist.fSeed[34],seed,sp);
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
        if (seed1->nHits() < 2) continue;

        for (int i2=i1+1; i2<nseeds; i2++) {

          DeltaSeed* seed2 = _data->deltaSeed(s,i2);
          if (seed2->nHits() < 2) continue;

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
      DeltaCandidate*      delta = _data->deltaCandidate(i);
      DeltaCandidatePar_t* dcp   = &_deltaPar[i];

      int   pdg_id = dcp->fMcPart->fPdgID;
      float mom    = dcp->fMcPart->Momentum();

      fillDeltaHistograms(_hist.fDelta[0],delta,dcp);

      if      (pdg_id > 2000)    fillDeltaHistograms(_hist.fDelta[1],delta,dcp);
      else if (pdg_id ==  11) {
        if (mom < 20)            fillDeltaHistograms(_hist.fDelta[2],delta,dcp);
      }
      else if (pdg_id      == -11)    fillDeltaHistograms(_hist.fDelta[6],delta,dcp);
      else if (abs(pdg_id) ==  13)    fillDeltaHistograms(_hist.fDelta[7],delta,dcp);
      else if (abs(pdg_id) == 211)    fillDeltaHistograms(_hist.fDelta[8],delta,dcp);

      if (delta->nHits() >= 5) {
        fillDeltaHistograms(_hist.fDelta[10],delta,dcp);

        if      (pdg_id      > 2000) fillDeltaHistograms(_hist.fDelta[11],delta,dcp);
        else if (pdg_id ==  11) {
          if (mom < 20)              fillDeltaHistograms(_hist.fDelta[12],delta,dcp);
        }
        else if (pdg_id      == -11) fillDeltaHistograms(_hist.fDelta[16],delta,dcp);
        else if (abs(pdg_id) ==  13) fillDeltaHistograms(_hist.fDelta[17],delta,dcp);
        else if (abs(pdg_id) == 211) fillDeltaHistograms(_hist.fDelta[18],delta,dcp);

        if (dcp->fNHitsCE > 0)       fillDeltaHistograms(_hist.fDelta[20],delta,dcp);
      }
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _listOfMcParticles.GetEntriesFast();

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = (McPart_t*) _listOfMcParticles.UncheckedAt(i);
      float mc_mom = mc->Momentum();

      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if (mc->fPdgID == PDGCode::e_minus) {
        fillMcHistograms(_hist.fMc[1],mc);

        if (mc->Time() > 550)             fillMcHistograms(_hist.fMc[2],mc);

        if (mc->fLastStation > mc->fFirstStation) {
          fillMcHistograms(_hist.fMc[3],mc);

          if (mc->NHits() >= 5)  {
            fillMcHistograms(_hist.fMc[4],mc);

            if (mc->fDelta != NULL) fillMcHistograms(_hist.fMc[5],mc);
            if (mc_mom     <    20) fillMcHistograms(_hist.fMc[6],mc);
          }
        }

        if ((mc_mom > 20) && (mc_mom < 80)) fillMcHistograms(_hist.fMc[100],mc);
        if (mc_mom > 80                   ) fillMcHistograms(_hist.fMc[101],mc);
      }
      else if (mc->fPdgID > 2000) {
//-----------------------------------------------------------------------------
// protons and such
// hits are the single-straw hits
//-----------------------------------------------------------------------------
        fillMcHistograms(_hist.fMc[200],mc);
        if (mc->NHits() > 10) fillMcHistograms(_hist.fMc[201],mc);
        if (mc->NHits() > 15) fillMcHistograms(_hist.fMc[202],mc);
      }
    }
//-----------------------------------------------------------------------------
// reconstructed protons
//-----------------------------------------------------------------------------
    int nprot = _data->nProtonCandidates();
    for (int i=0; i<nprot; i++) {
      ProtonCandidate*      pc  = _data->protonCandidate(i);
      ProtonCandidatePar_t* pcp = &_protonPar[i];

      fillProtonHistograms(_hist.fProton[0],pc,pcp);

      int nch = pc->nHitsTot();
      if (nch >= 10) fillProtonHistograms(_hist.fProton[1],pc,pcp);
      if (nch >= 15) fillProtonHistograms(_hist.fProton[2],pc,pcp);
    }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker
// the hope is that it is shorter than the list of all particles
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::initMcDiag() {
//-----------------------------------------------------------------------------
// memory cleanup after previous event
// assume that MC collections have been initialized
//-----------------------------------------------------------------------------
    int n = _listOfMcParticles.GetEntriesFast();
    for (int i=0; i<n; i++) {
      McPart_t* p = (McPart_t*) _listOfMcParticles.UncheckedAt(i);
      delete p;
    }

    _listOfMcParticles.Clear();
    _listOfMcPartHit.Clear  ();
    _listOfMcProtons.Clear  ();

    int   nch           = _data->chcol->size();
    if (nch <= 0) return -1;

    const ComboHit* ch0 = &_data->chcol->at(0);

    _listOfMcPartHit.Expand(nch);

    for (int ist=0; ist<kNStations; ist++) {
      for (int face=0; face<kNFaces; face++) {
        FaceZ_t* fz = &_data->fFaceData[ist][face];

        int nhits = fz->fHitData.size();
        for (int ih=0; ih<nhits; ih++) {
          HitData_t*      hd = &fz->fHitData[ih];
          const ComboHit* ch = hd->fHit;
          size_t ich         = ch-ch0;  // hit index in the collection
//-----------------------------------------------------------------------------
// associate combo hit with the MC particle hitting the first straw of a combo hit
// a combo hit could be made out of straw hits produced by different particles
// in neighboring strawa, so there could be an uncertainty - ignore that for now
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
            _listOfMcParticles.Add(mc);
            mc->fID       = _mcUtils->getID(sim);
            mc->fMotherID = _mcUtils->getMotherID(sim);
            mc->fPdgID    = _mcUtils->getPdgID(sim);
            mc->fStartMom = _mcUtils->getStartMom(sim);
            if (mc->fPdgID > 2000) {
                                        // add new proton
              _listOfMcProtons.Add(mc);
            }
          }
//-----------------------------------------------------------------------------
// list of hits produced by the particle
//-----------------------------------------------------------------------------
          mc->fListOfHits.push_back(hd);
//-----------------------------------------------------------------------------
// count N(hits) flagged as delta
//-----------------------------------------------------------------------------
          const StrawHitFlag* flag = &_data->outputChColl->at(ich).flag();
          if (flag->hasAnyProperty(StrawHitFlag::bkg)) mc->fNChFlaggedDelta += 1;

          int station        = ch->strawId().station();

          if (station < mc->fFirstStation) mc->fFirstStation = station;
          if (station > mc->fLastStation ) mc->fLastStation  = station;

          if (ch->correctedTime() < mc->fTime   ) mc->fTime    = ch->correctedTime();
          if (ch->correctedTime() < mc->fHitTMin) mc->fHitTMin = ch->correctedTime();
          if (ch->correctedTime() > mc->fHitTMax) mc->fHitTMax = ch->correctedTime();
//-----------------------------------------------------------------------------
// list of MC particles parallel to StrawHitCollection
//-----------------------------------------------------------------------------
          _listOfMcPartHit.AddAt(mc,ich);
        }
      }
    }

    _nMcProtons10 = 0;
    _nMcProtons15 = 0;
    _nMcProtons   = _listOfMcProtons.GetEntriesFast();

    for (int i=0; i<_nMcProtons; i++) {
      McPart_t* mcp = mcProton(i);
      if (mcp->NHits() >= 10) _nMcProtons10++;
      if (mcp->NHits() >= 15) _nMcProtons15++;
    }

    if (_data->debugLevel > 10) {
      int nmc = _listOfMcParticles.GetEntriesFast();
      printf(" N(MC particles with hits in the tracker: %5i\n",nmc);
      printf("    i     SimID        PdgID  NHits   Momentum  Time   FirstSt LastSt\n");
      for (int i=0; i<nmc; i++) {
        McPart_t* mc = (McPart_t*)_listOfMcParticles.At(i);
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
// initialize MC parameters of the reconstructed DeltaSeed's
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::initDeltaSeedMcParameters() {

    const ComboHit* hit0 = &_data->chcol->at(0);
//-----------------------------------------------------------------------------
// pre-calculate seed parameters
//-----------------------------------------------------------------------------
    for (int is=0; is<kNStations; is++) {
      int nseeds = _data->NSeeds(is);
      for (int i=0; i<nseeds; ++i) {
        DeltaSeed*      seed = _data->deltaSeed(is,i);
        DeltaSeedPar_t* sp   = &_seedPar[is].at(i);
//-----------------------------------------------------------------------------
// define MC pointers (SimParticle's) for the first two, "pre-seed", hits
//-----------------------------------------------------------------------------
        sp->fNHitsCE          = 0;
        sp->fPreSeedMcPart[0] = NULL;
        sp->fPreSeedMcPart[1] = NULL;

        int face0 = seed->SFace(0);
        if (seed->HitData(face0) != nullptr) {
          const HitData_t* hd1  = seed->HitData(face0);
          int loc1              = hd1->fHit-hit0;
          sp->fPreSeedMcPart[0] = mcPartHit(loc1);
        }

        int face1 = seed->SFace(1);
        if (face1 >= 0) {
          const HitData_t* hd2 = seed->HitData(face1);
          int loc2 = hd2->fHit-hit0;
          sp->fPreSeedMcPart[1] = mcPartHit(loc2);
        }
//-----------------------------------------------------------------------------
// define MC pointers (McPart_t's) for hits in all faces
// assume one ComboHit per face
//-----------------------------------------------------------------------------
        for (int face=0; face<kNFaces; face++) {
          const HitData_t* hd = seed->HitData(face);
          sp->fMcPart[face] = nullptr;
          if (hd == nullptr)                                          continue;
          int loc           = hd->fHit-hit0;         // ComboHit index in the collection
          McPart_t* mc      = mcPartHit(loc);
          sp->fMcPart[face] = mc;
//-----------------------------------------------------------------------------
// count CE hits
//-----------------------------------------------------------------------------
          if ((mc->fPdgID == PDGCode::e_minus) && (mc->fStartMom > 95) && (mc->fStartMom <110)) {
            sp->fNHitsCE += 1;
            if (seed->fDeltaIndex >= 0) {
              DeltaCandidatePar_t* dcp = &_deltaPar[seed->deltaIndex()];
              dcp->fNHitsCE           += 1;
            }
          }
        }
//-----------------------------------------------------------------------------
// count number of different simID's contributing
//-----------------------------------------------------------------------------
        sp->fNSim = 0;
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
            printSegment(seed,sp,0);
          }
        }
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
  int DeltaFinderDiag::initDeltaCandMcParameters() {
//-----------------------------------------------------------------------------
// proceed with the reconstructed delta candidates
// a) initialize additional delta parameters to be calculated
//-----------------------------------------------------------------------------
    int const    max_part(1000);
    McPart_t*    part [max_part];
    int          nhits[max_part];

    int ndelta = _data->nDeltaCandidates();
    _deltaPar.resize(ndelta);

    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidatePar_t* dcp = &_deltaPar[idelta];
      dcp->fMcPart   = nullptr;
      dcp->fNHitsMcP = 0;
      dcp->fNHitsCE  = 0;
    }

    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = _data->deltaCandidate(idelta);
//-----------------------------------------------------------------------------
// initialize residuals, whatever they are, to zero
//-----------------------------------------------------------------------------
      DeltaCandidatePar_t* dcp = &_deltaPar[idelta];

      float chi2_par, chi2_perp, xc, yc;

      xc = dc->Xc();
      yc = dc->Yc();

      _data->_finder->deltaChi2(dc,xc,yc,chi2_par,chi2_perp);

      dcp->fChi2ParN  = chi2_par /dc->nHits();
      dcp->fChi2PerpN = chi2_perp/dc->nHits();
      dcp->fChi2N     = dcp->fChi2ParN+dcp->fChi2PerpN;

      for (int is=0; is<kNStations; is++) dcp->dxy[is] = 0;

      int npart          = 0;
      for (int is=dc->fFirstStation; is<=dc->fLastStation; is++) {
        DeltaSeed* ds = dc->Seed(is);
        if (ds == 0)                                                  continue;

        DeltaSeedPar_t* sp = &_seedPar[is].at(ds->Index());

        for (int face=0; face<kNFaces; face++) {
          HitData_t* hd = (HitData_t*) ds->HitData(face);
          if (hd == nullptr)                                          continue;
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
              printf("* ERROR DeltaFinder::%s: npart >= max_part (%i)\n",__func__,max_part);
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
        dcp->fMcPart         = part [ibest];
        dcp->fNHitsMcP       = nhits[ibest];
        dcp->fMcPart->fDelta = dc;
      }
    }
    return 0;
  }


//-----------------------------------------------------------------------------
// proceed with the reconstructed proton candidates
// proton candidate is a time cluster - doesn't necessarily have to contain
// hits of only one particle
// however want to know !
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::initProtonCandMcParameters() {
    int const    max_part(1000);
    McPart_t*    part [max_part];
    int          nhits[max_part];

    const ComboHit* ch_hit_0 = &_data->chcol->at(0);

    _nRecoProtons = _data->nProtonCandidates();
    _protonPar.resize(_nRecoProtons);

    _nRecoProtons10 = 0;
    _nRecoProtons15 = 0;

    for (int i=0; i<_nRecoProtons; i++) {
      ProtonCandidatePar_t* pcp = &_protonPar[i];
      pcp->fMcPart   = nullptr;
      pcp->fNHitsMcP = 0;
      pcp->fNHitsCE  = 0;

      for (int is=0; is<kNStations; is++) {
        for (int face=0; face<kNFaces; face++) pcp->fHitMcPart[is][face].clear();
      }
    }

    for (int iprot=0; iprot<_nRecoProtons; iprot++) {
      ProtonCandidate*      pc  = _data->protonCandidate(iprot);
      ProtonCandidatePar_t* pcp = &_protonPar[iprot];
//-----------------------------------------------------------------------------
// two reconstructed proton counters - with 10+ and 15+ hits
//-----------------------------------------------------------------------------
      if (pc->nHitsTot() >= 10) _nRecoProtons10++;
      if (pc->nHitsTot() >= 15) _nRecoProtons15++;
//-----------------------------------------------------------------------------
// npart - the number of particles in the 'proton' time cluster
//-----------------------------------------------------------------------------
      int npart = 0;
      for (int is=pc->fFirstStation; is<=pc->fLastStation; is++) {
        int nh_station = pc->nHitsStation(is);
        if (nh_station == 0)                                          continue;
        for (int face=0; face<kNFaces; face++) {
//-----------------------------------------------------------------------------
// number of hits in this face
//-----------------------------------------------------------------------------
          int nh = pc->nHits(is,face);
          for (int ih=0; ih<nh; ih++) {
//-----------------------------------------------------------------------------
// assign 'proton' index to each hit in a proton time cluster
//-----------------------------------------------------------------------------
            HitData_t* hd = (HitData_t*) pc->hitData(is,face,ih);
            // already done : hd->fProtonIndex = iprot;
//-----------------------------------------------------------------------------
// try to identify a reconstructed proton with some MC particle
//-----------------------------------------------------------------------------
            int loc = hd->fHit-ch_hit_0;
            McPart_t* mc = mcPartHit(loc);  // list parallel to the list of CH

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
                printf("* DeltaFinder::%s ERROR: npart >= max_part (%i)\n",__func__,max_part);
              }
            }
          }
        }
      }
//-----------------------------------------------------------------------------
// look at the particles contributed to the cluster and determine
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
        pcp->fMcPart          = part [ibest];
        pcp->fNHitsMcP        = nhits[ibest];
        pcp->fMcPart->fProton = pc;
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
        _eventNumber = en;
        int rc = initMcDiag();
        if (rc < 0) return -1;
        associateMcTruth();
      }
    }

    if (_printComboHits) {
//-----------------------------------------------------------------------------
// print combo hits Z-ordered
//-----------------------------------------------------------------------------
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* DeltaFinderDiag::debug: all combo hit data\n");

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

    if (_printGoodComboHits) {
//-----------------------------------------------------------------------------
// print combo hits Z-ordered
//-----------------------------------------------------------------------------
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* DeltaFinderDiag::debug: good (non-delta and non-proton) combo hit data\n");

      const ComboHit* ch_hit_0 = &_data->chcol->at(0);

      for (int is=0; is<kNStations; is++) {
        printHitData(nullptr,1);
        for (int face=0; face<kNFaces; face++) {
          FaceZ_t* fz = &_data->fFaceData[is][face];
          int nhits = fz->fHitData.size();
          for (int ih=0; ih<nhits; ih++) {
            HitData_t* hd = &fz->fHitData[ih];
            int loc = hd->fHit-ch_hit_0;
            const StrawHitFlag* flag = &_data->outputChColl->at(loc).flag();

            if (  flag->hasAnyProperty(StrawHitFlag::bkg      )) continue;
            if (! flag->hasAnyProperty(StrawHitFlag::energysel)) continue;

            printHitData(hd,-1);
          }
        }
      }
    }

    if (_printDeltaSeeds != 0) {
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* DeltaFinderDiag::debug: print delta seeds\n");
      for (int st=0; st<kNStations; ++st) {
        int nseeds = _data->NSeeds(st);
        printf("* station: %2i N(seeds): %3i\n",st,nseeds);
        if (nseeds > 0) {

          printSegment(nullptr,nullptr,1);

          for (int iseed=0; iseed<nseeds; ++iseed) {
            DeltaSeed* seed = _data->deltaSeed(st,iseed);
            DeltaSeedPar_t* sp = &_seedPar[st].at(iseed);
            printSegment(seed,sp,-1);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// print reconstructed delta candidates
//-----------------------------------------------------------------------------
    if (_printDeltaCandidates != 0) printDeltaCandidates();
//-----------------------------------------------------------------------------
// print MC electrons
//-----------------------------------------------------------------------------
    if (_printElectrons) {
      int nmc = _listOfMcParticles.GetEntriesFast();

      for (int i=0; i<nmc; i++) {
        McPart_t* mc = (McPart_t*) _listOfMcParticles.At(i);

        if (((mc->fPdgID     == PDGCode::e_minus      ) or (mc->fPdgID     == PDGCode::e_plus)) &&
            (mc->Momentum() >  _printElectronsMinMom  ) &&
            (mc->Momentum() <  _printElectronsMaxMom  ) &&
            (mc->NHits()    >= _printElectronsMinNHits)    ) {

          float fr   = mc->fNHitsDelta/(mc->NHits()+1.e-3);
          float fbkg = mc->fNChFlaggedDelta/(mc->NHits()+1.e-3);

          if (fr < _printElectronsMaxFReco) {

            int delta_id(-1), nseeds(0), delta_nh(0), delta_nsh(0);
            if (mc->fDelta) {
              delta_id   = mc->fDelta->Index();
              nseeds     = mc->fDelta->NSeeds();
              delta_nh   = mc->fDelta->nHits();
              delta_nsh  = mc->fDelta->nStrawHits();
            }

            printf("* event: %4i electron.sim.id:%6i",_data->event->event(),mc->fID);
            printf(" mom:%7.3f time:%8.3f (deltaID: %3i nseg:nh:nsh:nf %2i:%2i:%2i)(nhits:nd:nf %3i:%2i:%2i) stations:%2i:%2i",
                   mc->Momentum(), mc->Time(),
                   delta_id, nseeds, delta_nh, delta_nsh,
                   mc->NHits(),
                   mc->fNHitsDelta,
                   mc->fNChFlaggedDelta,
                   mc->fFirstStation, mc->fLastStation);
            printf(" freco:fbkg %5.3f:%5.3f \n",fr,fbkg);

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
//-----------------------------------------------------------------------------
// print proton segments
//-----------------------------------------------------------------------------
    if (_printProtonSeeds != 0) {
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* DeltaFinderDiag::debug: print proton segments\n");
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      for (int st=0; st<kNStations; ++st) {
        int nseg = _data->nProtonSeeds(st);
        printf("* station: %2i N(proton segments): %3i\n",st,nseg);
        if (nseg > 0) {

          printSegment(nullptr,nullptr,1);

          for (int iseg=0; iseg<nseg; ++iseg) {
            DeltaSeed*      seg = _data->protonSeed(st,iseg);
            DeltaSeedPar_t* sp  = &_seedPar[st].at(seg->Index());
            printSegment(seg,sp,-1);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// print reconstructed proton candidates
// protons are more like time clusters
//-----------------------------------------------------------------------------
    if (_printProtonCandidates != 0) {
      const ComboHit* ch_hit_0 = &_data->chcol->at(0);
      int np = _data->nProtonCandidates();
      printf("* ----------------------------------------------------------------------------------------------------------------\n");
      printf("* [DeltaFinder::debug] N(proton candidates) = %5i\n",np);
      printf("* ----------------------------------------------------------------------------------------------------------------\n");

      if (np > 0) {
        for (int i=0; i<np; i++) {
          ProtonCandidate* pc       = _data->protonCandidate(i);
          ProtonCandidatePar_t* pcp = &_protonPar[i];
          int   pdg_id = -1;
          int   sim_id = -1;
          float mom    = -1.;
          if (pcp->fMcPart) {
            pdg_id = pcp->fMcPart->fPdgID;
            sim_id = pcp->fMcPart->fID;
            mom    = pcp->fMcPart->Momentum();
          }
          printf("* :pc:%05i itime:%3i",pc->index(),pc->timeIndex());
          printf("  nh:nCE %3i:%2i", pc->nHitsTot(),pcp->fNHitsCE);
          printf(" s1:s2 : %02i:%02i",pc->fFirstStation,pc->fLastStation);
          printf(" simID=%5i pdgID=%10i mom=%8.2f\n",sim_id,pdg_id,mom);
          printf("----------------------------------------------------------------------------------------------------------------\n");
          printf(" is  nh      Xc      Yc      Phi   Tmean     face=0       face=1       face=2      face=3\n");
          printf("----------------------------------------------------------------------------------------------------------------\n");

          for (int is=pc->fFirstStation;is<=pc->fLastStation; is++) {
            if (pc->nHitsStation(is) == 0) continue;

            printf("%3i  %3i ",is,pc->nHitsStation(is));
            printf("%7.1f %7.1f %7.2f %7.1f ",pc->xMean(is), pc->yMean(is), pc->phi(is), pc->t0(is));

            int nhits_face_max=0;
            for (int face=0; face<kNFaces; face++) {
              int nh = pc->nHits(is,face);
              if (nh > nhits_face_max) nhits_face_max = nh;
            }

            for (int i=0; i<nhits_face_max; i++) {
              if (i>0) printf("                                         ");
              for (int face=0; face<kNFaces; face++) {
                int nh = pc->nHits(is,face);
//-----------------------------------------------------------------------------
// form line to print
//-----------------------------------------------------------------------------
                if (nh <= i) printf("(           )");
                else {
                  // std::vector<McPart_t*>* hit_mcp = &pcp->fHitMcPart[is][face];
                  const HitData_t* hd = pc->hitData(is,face,i);
                  const ComboHit* hit = hd->fHit;
                  int loc = hit-ch_hit_0;
                  McPart_t* mcp = mcPartHit(loc);
                  int mcid      = mcp->fID;
                  printf("(%5i:%5i)",hit->strawId().asUint16(),mcid);
                }
              }
              printf(" \n");
            }
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// print MC protons and such .. PDG codes of all baryons are > 2000
//-----------------------------------------------------------------------------
    // constexpr PDGCode::enum_type kBaryon = PDGCode::enum_type(2000);

    if (_printMcProtons) {
      int nmc = _listOfMcParticles.GetEntriesFast();

      for (int i=0; i<nmc; i++) {
        McPart_t* mc = (McPart_t*) _listOfMcParticles.At(i);

        int pdg_id = (int) mc->fPdgID;
        if (pdg_id > 2000) {
          // &&
          //   (mc->Momentum() >  _printElectronsMinMom  ) &&
          //   (mc->Momentum() <  _printElectronsMaxMom  ) &&
          //   (mc->NHits()    >= _printElectronsMinNHits)    ) {

          float fr   = mc->fNHitsProton/(mc->NHits()+1.e-3);
          float fbkg = mc->fNChFlaggedProton/(mc->NHits()+1.e-3);

          if (fr < _printElectronsMaxFReco) {

            int proton_id(-1), proton_nseg(0), proton_nh(0), proton_nsh(0);
            if (mc->fProton) {
              proton_id    = mc->fProton->index();
              proton_nseg  = mc->fProton->nStationsWithHits();
              proton_nh    = mc->fProton->nHitsTot();
              proton_nsh   = mc->fProton->nStrawHitsTot();
            }

            printf("* event: %4i proton.sim.id:%6i",_data->event->event(),mc->fID);
            printf(" mom:%7.3f time:%8.3f (nhits:nd:nf %3i:%2i:%2i) (ProtonID:%5i nseg:nh:nsh:nf %2i:%2i:%2i) stations:%2i:%2i",
                   mc->Momentum(), mc->Time(),
                   mc->NHits(), mc->fNHitsProton, mc->fNChFlaggedProton,
                   proton_id, proton_nseg, proton_nh, proton_nsh,
                   mc->fFirstStation, mc->fLastStation);
            printf(" freco:fbkg %5.3f:%5.3f \n",fr,fbkg);

            if (_printProtonHits > 0) {
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
// this function also resizes the array of delta candidate parameters
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::precalculateRecoParameters() {

    for (int is=0; is<kNStations; is++) {
      int nseeds = _data->NSeeds(is);
      _seedPar[is].resize(nseeds);
    }

    for (int s=kNStations-1; s>=0; s--) {
      int nseeds = _data->NSeeds(s);
      for(int se=0; se<nseeds; ++se) {
        DeltaSeed*      seed = _data->deltaSeed(s,se);
        DeltaSeedPar_t* sp   = &_seedPar[s].at(se);

        float xc = seed->Xc();
        float yc = seed->Yc();

        sp->fDt12     = -1.e6;
        sp->fDtCorr12 = -1.e6;
        sp->fDtDelta  = -1.e6;

        if (seed->SFace(1) >= 0) {
          sp->fDt12     = seed->Dt12();
          sp->fDtCorr12 = seed->DtCorr12();
        }

        for (int face=0; face<kNFaces; face++) {
          HitData_t* hd = seed->HitData(face);
          if (hd) {
//-----------------------------------------------------------------------------
// use DeltaFinderAlg to calculate chi2
//-----------------------------------------------------------------------------
            _data->_finder->hitChi2(hd,xc,yc,sp->fChi2HPar[face],sp->fChi2HPerp[face]);
            sp->fChi2H[face] = sp->fChi2HPar[face] + sp->fChi2HPerp[face];
          }
          else {
            sp->fChi2H    [face] = -1;
            sp->fChi2HPar [face] = -1;
            sp->fChi2HPerp[face] = -1;
          }
        }
        if (seed->Chi2ParN() > _data->_finder->_maxChi2Par) {
          if (_data->_finder->printErrors()) {
            printf("* ERROR .. should never happen! seed chi2par/N = %10.3f\n",seed->Chi2ParN());
            printSegment(seed,sp,0);
          }
        }
      }
    }

    int ndelta = _data->nDeltaCandidates();
    _deltaPar.resize(ndelta);

    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate*      dc  = _data->deltaCandidate(idelta);
      DeltaCandidatePar_t* dcp = &_deltaPar[idelta];
//-----------------------------------------------------------------------------
// initialize residuals, whatever they are, to zero
//-----------------------------------------------------------------------------
      for (int is=0; is<kNStations; is++) dcp->dxy[is] = 0;

      for (int is=dc->fFirstStation; is<=dc->fLastStation; is++) {
        DeltaSeed* ds = dc->Seed(is);
        if (ds == 0)                                                  continue;
//-----------------------------------------------------------------------------
// to calculate the seed residuals wrt the delta, exclude the seed
// so 'dc1' = 'dc' - 'ds'
//-----------------------------------------------------------------------------
        DeltaSeedPar_t* sp = &_seedPar[is].at(ds->Index());

        DeltaCandidate dc1;
        dc1.removeSeed(dc,is);
//-----------------------------------------------------------------------------
// calculate the seed residuals and chi2 wrt dc1
//-----------------------------------------------------------------------------
        float    chi2_par, chi2_perp;
        _data->_finder->seedChi2(ds,dc1.Xc(),dc1.Yc(),chi2_par,chi2_perp);

        sp->fDeltaChi2Par  = chi2_par /ds->nHits();
        sp->fDeltaChi2Perp = chi2_perp/ds->nHits();
        sp->fDeltaChi2     = (chi2_par+chi2_perp)/ds->nHits();
        sp->fDtDelta       = ds->TMean() - dc1.T0(is);
      }
    }
    return 0;
  }
//-----------------------------------------------------------------------------
// print flag stored in the ComboHit output hit flag collection
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printHitData(const HitData_t* Hd, int PrintBanner) {

    if (PrintBanner >= 0) {
      printf("#----------------------------------------------------------------------------------------------");
      printf("---------------------------------------------------------------------------------------------\n");
      printf("#   I  SHID    flag   nsh St:F:Pl:Pn:Str   Time     TCorr     dt     eDep     wdist     wres   ");
      printf("     PDG simID    mom_PDG mom_ID   pStart      p       pz       X        Y         Z  DeltaID ProtID\n");
      printf("#----------------------------------------------------------------------------------------------");
      printf("---------------------------------------------------------------------------------------------\n");
      if (PrintBanner > 0) return;
    }

    const ComboHit* ch0 = &_data->chcol->at(0);
    const ComboHit* ch  = Hd->fHit;
    int loc             = ch-ch0;

    const StrawHitFlag* flag = &(*_data->outputChColl)[loc].flag();

    // int radselOK        = (! flag->hasAnyProperty(StrawHitFlag::radsel   ));
    // int edepOK          = (! flag->hasAnyProperty(StrawHitFlag::energysel));

    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999), mom_id(-9999), mom_pdg_id(-9999);
    float              mc_mom(-9999.);
    const XYZVectorF*  p(nullptr);

    if (_mcDiag) {
      sim        = _mcUtils->getSimParticle(_data->event,ch->index(0));
      p          = _mcUtils->getMom(_data->event,ch->index(0));
      pdg_id     = _mcUtils->getPdgID(sim);
      sim_id     = _mcUtils->getID(sim);
      mc_mom     = _mcUtils->getStartMom(sim);
      mom_id     = _mcUtils->getMotherID(sim);
      mom_pdg_id = _mcUtils->getMotherPdgID(sim);
    }

    printf("%5i",loc);
    printf(" %5i 0x%08x %2i" ,ch->strawId().asUint16(),*((int*) flag),ch->nStrawHits());

    printf(" %2i:%i:%02i %2i %2i %8.2f %8.2f %7.3f %8.5f %8.3f %8.3f %10i %5i %10i %5i %8.3f %8.3f %8.3f %8.3f %8.3f %9.3f %5i %6i\n",
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
           p->r(),
           p->z(),
           ch->pos().x(),
           ch->pos().y(),
           ch->pos().z(),
           Hd->DeltaIndex(),
           Hd->ProtonIndex()
           );
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

    const StrawHitFlag* shf = &Ch->flag(); // _data->outputChfColl->at(loc);

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
  void DeltaFinderDiag::printDeltaCandidates() {
    int nd = _data->nDeltaCandidates();
    printf("* ----------------------------------------------------------------------------------------------------------------\n");
    printf("* [DeltaFinder::debug] N(delta candidates) = %5i\n",nd);
    printf("* ----------------------------------------------------------------------------------------------------------------\n");

    if (nd > 0) {
      for (int i=0; i<nd; i++) {
        DeltaCandidate* dc       = _data->deltaCandidate(i);
        DeltaCandidatePar_t* dcp = &_deltaPar[i];
        int pdg_id = -1;
        int sim_id = -1;
        if (dcp->fMcPart) {
          pdg_id = dcp->fMcPart->fPdgID;
          sim_id = dcp->fMcPart->fID;
        }

        printf("* :dc:%05i",dc->Index());
        printf("  nh n(CE) ns s1  s2     X        Y        Z       R     phi chi21 chi22  htmin  htmax   t0   PdgID N(MC hits)");
        printf(" simID=%5i pdgID=%10i \n",sim_id,pdg_id);
        printf("----------------------------------------------------------------------------------------------------------------\n");
        printf("            %3i  %3i %3i",dc->fNHits,dcp->fNHitsCE,dc->fNSeeds);
        printf(" %2i  %2i %7.2f %7.2f %9.2f",dc->fFirstStation,dc->fLastStation,
               dc->CofM.x(),dc->CofM.y(),dc->CofM.z());
        printf("               %39i %5i",pdg_id,dcp->fNHitsMcP);
        printf("\n");
        printf("----------------------------------------------------------------------------------------------------------------\n");

        for (int is=dc->fFirstStation;is<=dc->fLastStation; is++) {
          DeltaSeed* ds = dc->Seed(is);
          if (ds != NULL) {

            DeltaSeedPar_t* sp = &_seedPar[is].at(ds->Index());

            int face0 = ds->SFace(0);
            int face1 = ds->SFace(1);

            const HitData_t* hd0 = ds->HitData(face0);
            const HitData_t* hd1 = (face1 >= 0) ? ds->HitData(face1) : nullptr;

            float x = ds->CofM.x();
            float y = ds->CofM.y();
            float z = ds->CofM.z();

            if (ds->nHits() == 1) {
              x = hd0->fHit->pos().x();
              y = hd0->fHit->pos().y();
              z = hd0->fHit->pos().z();
            }

            float r   = sqrt(x*x+y*y);
            float phi = atan2(y,x);

            printf("            %3i  %3i    %3i:%03i",ds->fNHits,sp->fNHitsCE,is,ds->Index());
            printf(" %7.2f %7.2f %9.2f %7.2f %5.2f",x,y,z,r,phi);
            float chi22 = (hd1) ? hd1->fChi2Min : -1;
            printf(" %5.1f %5.1f",hd0->fChi2Min, chi22);
            printf(" %6.1f %6.1f",ds->MinHitTime(),ds->MaxHitTime());
            printf(" %6.1f ",dc->T0(is));

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
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printSegment(DeltaSeed* seed, DeltaSeedPar_t* Sp, int PrintBanner) {

    if (PrintBanner >= 0) {
      printf("* -----------------------------------------------------------------------------------------");
      printf("-------------------------------------------------------------------------------\n");
      printf("* st seed  good type delta proton SHID: MCID   SHID: MCID   SHID: MCID   SHID: MCID");
      printf("   chi21  chi22 chi2all/N chipar/N ch2prp/N mintime  maxtime  <edep>      ");
      printf("X        Y         Z   nch nsh\n");
      printf("* -------------------------------------------------------------------------------------------");
      printf("-------------------------------------------------------------------------------\n");
    }

    if (PrintBanner > 0) return;

    printf("* %2i  %03i %5i   %02i",seed->Station(),seed->Index(),seed->fGood,seed->fType);
//-----------------------------------------------------------------------------
// it is important to print seed->fDeltaIndex, not seed->deltaIndex(),
// as the latter performs the truncation: fDeltaIndex % 10000
//-----------------------------------------------------------------------------
    printf(" %5i %5i",seed->deltaIndex(),seed->protonIndex());
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
    for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd = seed->HitData(face);
      if (hd == nullptr) printf("(%5i:%5i)",-1,-1);
      else {
        const ComboHit* hit = hd->fHit;
        McPart_t* mcp = Sp->fMcPart[face];
        if (mcp) printf("(%5i:%5i)",hit->strawId().asUint16(),mcp->fID);
        else     printf("(%5i:XXXXX)",hit->strawId().asUint16());

      }
    }

    printf(" %6.2f %6.2f %8.2f %8.2f %8.2f",seed->fChi21,seed->fChi22, seed->Chi2TotN(),seed->Chi2ParN(),seed->Chi2PerpN());
    printf(" %8.1f %8.1f %8.5f",seed->MinHitTime(),seed->MaxHitTime(),seed->EDep());
    printf(" %8.3f %8.3f %9.3f",seed->CofM.x(),seed->CofM.y(),seed->CofM.z());
    //    printf("%4i",seed->fNFacesWithHits);
    printf("%4i",seed->nHits());
    printf("%4i",seed->nStrawHits());
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
