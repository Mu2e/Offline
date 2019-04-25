#ifndef __CalPatRec_DeltaFinderDiag_hh__
#define __CalPatRec_DeltaFinderDiag_hh__

#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "CalPatRec/inc/DeltaFinder_types.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Mu2eUtilities/inc/McUtilsToolBase.hh"

using namespace std;

namespace mu2e {
  using namespace DeltaFinderTypes;

  class SimParticle;
  class StrawDigiMCCollection;
  
  class DeltaFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets =  10,
      kNSeedHistSets  =  10,
      kNDeltaHistSets =  10,
      kNMcHistSets    = 200
    };

    struct SeedHist_t {
      TH1F*  fChi2N;
      TH1F*  fChi2Tot;
      TH1F*  fHitChi2Min;			// chi2 of the first two hits along the wire
      TH1F*  fChi2Neighbour;
      TH1F*  fChi2Radial;
      TH1F*  fNFacesWithHits;
      TH1F*  fNHitsPerFace;
      TH1F*  fNHitsPerSeed;
      TH1F*  fSeedRadius;
      TH1F*  fSeedMomentum;
      TH2F*  fSeedSize;
      TH1F*  fNHits;
    };

    struct DeltaHist_t {
      TH1F*  fNHits;
      TH1F*  fNSeeds;
      TH1F*  fMcMom;
      TH1F*  fPDGCode;
      TH1F*  fDxy;
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
      TH1F*  fNSecondHits;
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
      DeltaHist_t* fDelta[kNDeltaHistSets];
      McHist_t*    fMc   [kNMcHistSets   ];
    };
  protected:

    bool                                  _mcDiag;
    art::InputTag                         _stepPointMcCollTag;
    int                                   _printOTracker;
    int                                   _printElectrons;
    int                                   _printElectronsHits;
    int                                   _printElectronsMinNHits;
    float                                 _printElectronsMaxFReco;
    float                                 _printElectronsMinMom;
    float                                 _printElectronsMaxMom;
    int                                   _printDeltaSeeds;
    int                                   _printDeltaCandidates;
    int                                   _printShcol;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    int                                   _eventNumber;
    //    const StrawDigiMCCollection*          _listOfMcStrawHits;
    int                                   _nDeltaHitsTot;
    int                                   _nDeltaHitsReco;
    
    std::vector<McPart_t*>                _list_of_mc_particles; // list_of_particles with hits in the tracker
    std::vector<McPart_t*>                _list_of_mc_part_hit ; // for each StrawHit, pointer to its McPart 

    Data_t*                               _data;                 // diag data, passed from the caller, cached

    Hist_t                                _hist;

  public:
    
    DeltaFinderDiag(const fhicl::ParameterSet& PSet);
    ~DeltaFinderDiag();

  private:

    void        bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir);
    void        bookSeedHistograms (SeedHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookDeltaHistograms(DeltaHist_t* Hist, art::TFileDirectory* Dir);
    void        bookMcHistograms   (McHist_t*    Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms(EventHist_t* Hist);
    void        fillSeedHistograms (SeedHist_t*  Hist, DeltaSeed*      Seed );
    void        fillDeltaHistograms(DeltaHist_t* Hist, DeltaCandidate* Delta);
    void        fillMcHistograms   (McHist_t*    Hist, McPart_t*       Mc   );

    McPart_t*   findParticle (const SimParticle* Sim);
    int         InitMcDiag      ();
    int         associateMcTruth();
    
    void        printHitData(const HitData_t* Hd, int Index);
    void        printOTracker();

    void        printComboHit(const ComboHit* Sh, int Index);
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
    _stepPointMcCollTag    (PSet.get<string>       ("stepPointMcCollTag"           )),
    _printOTracker         (PSet.get<int>          ("printOTracker"                )),
    _printElectrons        (PSet.get<int>          ("printElectrons"               )),
    _printElectronsHits    (PSet.get<int>          ("printElectronsHits"           )),
    _printElectronsMinNHits(PSet.get<int>          ("printElectronsMinNHits"       )),
    _printElectronsMaxFReco(PSet.get<float>        ("printElectronsMaxFReco"       )),
    _printElectronsMinMom  (PSet.get<float>        ("printElectronsMinMom"         )),
    _printElectronsMaxMom  (PSet.get<float>        ("printElectronsMaxMom"         )),
    _printDeltaSeeds       (PSet.get<int>          ("printDeltaSeeds"              )),
    _printDeltaCandidates  (PSet.get<int>          ("printDeltaCandidates"         )),
    _printShcol            (PSet.get<int>          ("printShcol"                   ))
  {
    printf(" DeltaFinderDiag::DeltaFinderDiag : HOORAY! \n");
    //    _firstCall   =  1;
    //    _timeOffsets = NULL;

    if (_mcDiag != 0) _mcUtils = art::make_tool<McUtilsToolBase>(PSet.get<fhicl::ParameterSet>("mcUtils"));
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }
  
//-----------------------------------------------------------------------------
  DeltaFinderDiag::~DeltaFinderDiag() {
  }
  
//-----------------------------------------------------------------------------
  McPart_t* DeltaFinderDiag::findParticle(const SimParticle* Sim) {
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
  void DeltaFinderDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number", 100, 0., 100000.);
    Hist->fNSecondHits     = Dir->make<TH1F>("nhit2" , "N(second hits)", 100, 0., 100.);
    Hist->fNSeeds          = Dir->make<TH1F>("nseeds", "N(seeds)"   , 200, 0., 2000.);
    Hist->fNSeedsVsStation = Dir->make<TH2F>("ns_vs_st", "N(seeds) vs station", 20, 0., 20.,100,0,100);

    Hist->fNMc             = Dir->make<TH1F>("nmc"       , "N(MC particles)", 100, 0., 1000.);
    Hist->fPDGCode         = Dir->make<TH1F>("pdg_code"  , "PDG Code"       , 100, 0., 100.);
    Hist->fNMcHits         = Dir->make<TH1F>("n_mc_hits" , "N(MC hits)"     , 250, 0., 500.);
    Hist->fNDelta          = Dir->make<TH1F>("n_delta"   , "N(reco deltas)" , 250, 0., 250.);
    Hist->fNDeltaHitsT     = Dir->make<TH1F>("n_delta_ht", "N(delta hits T)", 500, 0., 5000.);
    Hist->fNDeltaHitsR     = Dir->make<TH1F>("n_delta_hr", "N(delta hits R)", 500, 0., 5000.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookSeedHistograms(SeedHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fChi2N           = Dir->make<TH1F>("chi2n"       , "Chi2/N"            , 1000, 0., 100.);
    Hist->fChi2Tot         = Dir->make<TH1F>("chi2t"       , "Chi (all 2 hits)"  , 1000, 0., 100.);
    Hist->fHitChi2Min      = Dir->make<TH1F>("hitchi2min"  , "Hit Chi (min)"     , 1000, 0.,  50.);
    Hist->fChi2Neighbour   = Dir->make<TH1F>("chi2_nb"     , "Chi2 neighbour"    , 1000, 0.,  50.);
    Hist->fChi2Radial      = Dir->make<TH1F>("chi2_r"      , "Chi2 radial"       , 1000, 0.,  50.);
    Hist->fNFacesWithHits  = Dir->make<TH1F>("nfaces_wh"   , "Number of faces with hits", 5, 0., 5.);
    Hist->fNHitsPerFace    = Dir->make<TH1F>("nhits_face"  , "Number of hits per face", 20, 0., 20.);
    Hist->fNHitsPerSeed    = Dir->make<TH1F>("nhits_seed"  , "Number of hits per seed", 40, 0., 40.);
    Hist->fSeedRadius      = Dir->make<TH1F>("rad"         , "Seed radius", 500, 0., 1000.);
    Hist->fSeedMomentum    = Dir->make<TH1F>("mom"         , "Seed momentum", 400, 0., 400.);
    Hist->fSeedSize        = Dir->make<TH2F>("seed_size"   , "Seed (nh2+1):(nh1+1)", 20, 0., 20.,20,0,20);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookDeltaHistograms(DeltaHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fNSeeds  = Dir->make<TH1F>("nseeds", "N(seeds)" ,  20, 0.,  20.);
    Hist->fNHits   = Dir->make<TH1F>("nhits" , "N(hits)"  , 200, 0., 200.);
    Hist->fPDGCode = Dir->make<TH1F>("pdg"   , "PDG code" , 200, 0., 200.);
    Hist->fMcMom   = Dir->make<TH1F>("mc_mom", "N(hits)"  , 200, 0., 200.);
    Hist->fDxy     = Dir->make<TH1F>("dxy"   , "Delta Dxy", 500, 0., 200.);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::bookMcHistograms(McHist_t* Hist, art::TFileDirectory* Dir) {

    Hist->fPDGCode    = Dir->make<TH1F>("pdg"  , "PDG code"        , 500, -250., 250.);
    Hist->fMom        = Dir->make<TH1F>("mom"  , "momentum"        , 200, 0., 200.);
    Hist->fNHits      = Dir->make<TH1F>("nhits", "N(hits)"         , 200, 0., 200.);
    Hist->fNHitsDelta = Dir->make<TH1F>("nhitsr", "N(hits reco)"   , 200, 0., 200.);
    Hist->fFractReco  = Dir->make<TH1F>("fractr", "NR/N"           , 100, 0.,   1.);
    Hist->fMaxSeg     = Dir->make<TH1F>("max_seg", "Max N Segments",  20, 0.,  20.);
    Hist->fHitDt      = Dir->make<TH1F>("hit_dt" , "Hit TMax-TMin" , 100, 0., 200.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
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

    book_event_histset[ 0] = 1;		// all events

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

    book_seed_histset[ 0] = 1;		// all seeds
    book_seed_histset[ 1] = 0;          // *** unused
    book_seed_histset[ 2] = 1;          // truth 
    book_seed_histset[ 3] = 0;          // *** unused 
    book_seed_histset[ 4] = 1;          // truth and P < 30 MeV/c
    book_seed_histset[ 5] = 1;          // truth and P > 30 MeV/c
    book_seed_histset[ 6] = 1;          // truth = 0

    for (int i=0; i<kNSeedHistSets; i++) {
      if (book_seed_histset[i] != 0) {
	sprintf(folder_name,"seed_%i",i);
	art::TFileDirectory dir = Tfs->mkdir(folder_name);
	
	_hist.fSeed[i] = new SeedHist_t;
	bookSeedHistograms(_hist.fSeed[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book delta histograms
//-----------------------------------------------------------------------------
    int book_delta_histset[kNDeltaHistSets];
    for (int i=0; i<kNDeltaHistSets; i++) book_delta_histset[i] = 0;

    book_delta_histset[ 0] = 1;		// all  deltas
    book_delta_histset[ 1] = 1;		// long deltas

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

    book_mc_histset[  0] = 1;		// all particles
    book_mc_histset[  1] = 1;		// electrons
    book_mc_histset[  2] = 1;		// electrons fTime > 550
    book_mc_histset[  3] = 1;		// electrons fTime > 550 with last>first
    book_mc_histset[  4] = 1;		// electrons fTime > 550 with last>first and 4+ hits
    book_mc_histset[  5] = 1;		// electrons fTime > 550 with last>first, 5+ hits, and reco delta
    book_mc_histset[  6] = 1;		// electrons fTime > 550 with last>first, 5+ hits, and p < 20

    book_mc_histset[100] = 1;		// electrons fTime > 550 and 20 < p < 80 MeV/c
    book_mc_histset[101] = 1;		// electrons fTime > 550 and p > 80 MeV/c

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
  void  DeltaFinderDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();
    
    Hist->fEventNumber->Fill(event_number);
    Hist->fNSeeds->Fill(_data->nseeds);

    for (int is=0; is<kNStations; ++is) {

      Hist->fNSeedsVsStation->Fill(is,_data->nseeds_per_station[is]);

      for (int f=0; f<kNFaces-1; ++f) { // loop over faces
	for (int p=0; p<kNPanelsPerFace; ++p) { // loop over panels
	  PanelZ_t* panelz = &_data->oTracker[is][f][p];

	  // for (int l=0; l<2; ++l) { //loop over layers
	  int nhits =  panelz->fHitData.size();
	  for (int i=0; i<nhits; i++) {
	    int counter  = panelz->fHitData.at(i).fNSecondHits;
	    Hist->fNSecondHits->Fill(counter);
	  }
	  // }
	}
      }
    }

    int ndelta = _data->deltaCandidateHolder.size();
    Hist->fNDelta->Fill(ndelta);
//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      int n_mc_hits = mc->fListOfHits.size();

      Hist->fNMcHits->Fill(n_mc_hits);
      Hist->fPDGCode->Fill(mc->fPdgID);
    }

    Hist->fNMc->Fill(nmc);
    Hist->fNDeltaHitsT->Fill(_nDeltaHitsTot);
    Hist->fNDeltaHitsR->Fill(_nDeltaHitsReco);
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderDiag::fillSeedHistograms(SeedHist_t* Hist, DeltaSeed* Seed) {

    DeltaFinderTypes::Intersection_t res;

    Hist->fChi2N->Fill(Seed->Chi2N());
    Hist->fChi2Tot->Fill(Seed->Chi2Tot());

    if (Seed->fType != 0) {
      int face1 = Seed->fType/10;
      int face2 = Seed->fType % 10;

      int nh1 = Seed->hitlist[face1].size();
      int nh2 = Seed->hitlist[face2].size();

      const HitData_t* hd1 = Seed->hitlist[face1].at(0);
      const HitData_t* hd2 = Seed->hitlist[face2].at(0);

      DeltaFinderTypes::findIntersection(hd1,hd2,&res);

      double chi1 = res.wd1/hd1->fSigW;
      double chi2 = res.wd2/hd2->fSigW;

      Hist->fHitChi2Min->Fill(chi1*chi1);
      Hist->fHitChi2Min->Fill(chi2*chi2);
      Hist->fSeedSize->Fill(nh1,nh2);

      for (int i=1; i<nh2; i++) {
	const HitData_t* hd  = Seed->hitlist[face2][i];
	DeltaFinderTypes::findIntersection(hd,hd1,&res);
	float chi = res.wd1/hd->fSigW;
	Hist->fChi2Neighbour->Fill(chi*chi);
      }

      for (int i=1; i<nh1; i++) {
	const HitData_t* hd  = Seed->hitlist[face1][i];
	DeltaFinderTypes::findIntersection(hd,hd2,&res);
	float chi = res.wd1/hd->fSigW;
	Hist->fChi2Neighbour->Fill(chi*chi);
      }

      double _sigmaR = 10;

      for (int face=0; face<kNFaces; face++) {
	if ((face == face1) || (face == face2)) continue;
	int nh = Seed->hitlist[face].size();
	for (int ih=0; ih<nh; ih++) {
	  const HitData_t* hd  = Seed->hitlist[face][ih];
//-----------------------------------------------------------------------------
// reproduce DeltaFinder algorithm
//-----------------------------------------------------------------------------
	  // const StrawHitPosition* shp  = hd->fPos;
	  CLHEP::Hep3Vector       dxyz = hd->fHit->posCLHEP()-Seed->CofM;//shp->posCLHEP()-Seed->CofM; // distance from hit to preseed
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
	  const CLHEP::Hep3Vector& wdir = hd->fHit->wdirCLHEP();//hd->fStraw->getDirection();
	  CLHEP::Hep3Vector d_par    = (dxyz.dot(wdir))/(wdir.dot(wdir))*wdir; 
	  CLHEP::Hep3Vector d_perp_z = dxyz-d_par;
	  float  d_perp              = d_perp_z.perp();
	  double sigw                = hd->fSigW;
	  float  chi2_par            = (d_par.mag()/sigw)*(d_par.mag()/sigw);
	  float  chi2_perp           = (d_perp/_sigmaR)*(d_perp/_sigmaR);
	  float  chi2r               = chi2_par + chi2_perp;
	  Hist->fChi2Radial->Fill(chi2r);
	}
      }
    }

    Hist->fNFacesWithHits->Fill(Seed->fNFacesWithHits);
    Hist->fSeedRadius->Fill    (Seed->CofM.perp());

    double mom (-1.);
    if (Seed->fPreSeedMcPart[0]) mom = Seed->fPreSeedMcPart[0]->Momentum();

    Hist->fSeedMomentum->Fill(mom);

    Hist->fNHitsPerSeed->Fill(Seed->NHitsTot());  
  }


//-----------------------------------------------------------------------------
  void DeltaFinderDiag::fillDeltaHistograms(DeltaHist_t* Hist, DeltaCandidate* Delta) {
    int n_seeds = Delta->n_seeds;
    Hist->fNSeeds->Fill(n_seeds);
    Hist->fNHits->Fill(Delta->fNHits);

    float mom(-1), pdg_code(-1.e6);
    if (Delta->fMcPart) {
      mom      = Delta->fMcPart->Momentum();
      pdg_code = Delta->fMcPart->fPdgID;
    }

    Hist->fPDGCode->Fill(pdg_code);
    Hist->fMcMom->Fill(mom);

    for(int is=0; is<kNStations; ++is) {
      Hist->fDxy->Fill(Delta->dxy[is]);
    }
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
	  _eventNumber       = en;
	  //	_listOfMcStrawHits = _mcUtils->getListOfMcStrawHits(_data->event, _stepPointMcCollTag);
	InitMcDiag();
	associateMcTruth();
      }
    }
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// per-seed histograms
//-----------------------------------------------------------------------------
    for (int s=0; s<kNStations; ++s) {
      int seedsize = _data->seedHolder[s].size();
      for(int se=0; se<seedsize; ++se) {
	DeltaSeed* seed = _data->seedHolder[s].at(se);
	
	fillSeedHistograms(_hist.fSeed[0],seed);

	if (seed->MCTruth()) {
//-----------------------------------------------------------------------------
// real pre-seed - made out of hits produced by the same particle
//-----------------------------------------------------------------------------
	  fillSeedHistograms(_hist.fSeed[2],seed);

	  float mom = seed->fPreSeedMcPart[0]->Momentum();
	  if (mom < 20) fillSeedHistograms(_hist.fSeed[4],seed);
	  else          fillSeedHistograms(_hist.fSeed[5],seed);
	}
	else {
//-----------------------------------------------------------------------------
// fake pre-seed - made out of hits produced by two different particles
//-----------------------------------------------------------------------------
	  fillSeedHistograms(_hist.fSeed[6],seed);
	}
      }
    }
//-----------------------------------------------------------------------------
// fill delta histograms
//-----------------------------------------------------------------------------
    int ndelta = _data->deltaCandidateHolder.size();
    for(int i=0; i<ndelta; i++) {
      DeltaCandidate* delta = &_data->deltaCandidateHolder.at(i);

      fillDeltaHistograms(_hist.fDelta[0],delta);
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);

      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if (mc->fPdgID == 11) {
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
// create a list of MC particles with hits in the tracker- hope, it is shorter
// than the list of all particles
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::InitMcDiag() {
//-----------------------------------------------------------------------------
// memory cleanup after previous event
// assume that MC collections have been initialized 
//-----------------------------------------------------------------------------
    int n = _list_of_mc_particles.size();
    for (int i=0; i<n; i++) {
      McPart_t* p = _list_of_mc_particles.at(i);
      delete p;
    }
    
    _list_of_mc_particles.clear();
    _list_of_mc_part_hit.clear();
    
    const ComboHit* sh0 = &_data->chcol->at(0);
    int   nsh           = _data->chcol->size();

    _list_of_mc_part_hit.resize(nsh);
    
    for (int ist=0; ist<kNStations; ist++) {
      for (int face=0; face<kNFaces; face++) {
	for (int ip=0; ip<kNPanelsPerFace; ip++) {
	  PanelZ_t* panelz = &_data->oTracker[ist][face][ip];
	  
	  // for (int il=0; il<2; il++) {
	  int nhits = panelz->fHitData.size();
	  for (int ih=0; ih<nhits; ih++) {
	    HitData_t* hd = &panelz->fHitData[ih];

	    const ComboHit* sh           = hd->fHit;
	    size_t ish                   = sh-sh0;
	    // get the StrawDigi indices associated with this ComboHit
	    std::vector<StrawDigiIndex> shids;
	    _data->chcol->fillStrawDigiIndices(*(_data->event),ish,shids);
	    const mu2e::SimParticle* sim = _mcUtils->getSimParticle(_data->event,shids[0]);//ish);
//-----------------------------------------------------------------------------
// search if this particle has already been registered
//-----------------------------------------------------------------------------
	    McPart_t* mc = findParticle(sim);

	    if (mc == NULL) {
	      // add new particle
	      mc = new McPart_t(sim);
	      _list_of_mc_particles.push_back(mc);
	      mc->fID       = _mcUtils->getID(sim);
	      mc->fPdgID    = _mcUtils->getPdgID(sim);
	      mc->fStartMom = _mcUtils->getStartMom(sim);
	    }
	    //-----------------------------------------------------------------------------
	    // list of hits produced by the particle
	    //-----------------------------------------------------------------------------
	    mc->fListOfHits.push_back(hd);
	      
	      // StrawId   shid  = sh->strawId();
	      // const Straw& straw = _data->tracker->getStraw(shid);
	      // int station        = straw.id().getStation();
	    int station        = sh->strawId().station();

	    if (station < mc->fFirstStation) mc->fFirstStation = station;
	    if (station > mc->fLastStation ) mc->fLastStation  = station;
	      
	    if (sh->time() < mc->fTime   ) mc->fTime    = sh->time();
	    if (sh->time() < mc->fHitTMin) mc->fHitTMin = sh->time();
	    if (sh->time() > mc->fHitTMax) mc->fHitTMax = sh->time();
	    //-----------------------------------------------------------------------------
	    // list of MC particles parallel to StrawHitCollection
	    //-----------------------------------------------------------------------------
	    _list_of_mc_part_hit[ish] = mc;
	  }
	  // }
	}
      }
    }

    if (_data->debugLevel > 10) {
      int nmc = _list_of_mc_particles.size();
      printf(" N(MC particles with hits in the tracker: %5i\n",nmc);
      printf("    i     SimID        PdgID  NHits   Momentum  Time   FirstSt LastSt\n");
      for (int i=0; i<nmc; i++) {
	McPart_t* mc = _list_of_mc_particles.at(i);
	printf(" %4i  %10i %10i  %5i %10.3f %8.1f %5i %5i\n",
	       i,mc->fID,mc->fPdgID,
	       mc->NHits(),
	       mc->Momentum(),
	       mc->Time(),
	       mc->fFirstStation,mc->fLastStation);
      }
    }
    return 0;
  }

//-----------------------------------------------------------------------------
// for each DeltaSeed, create a list of SimParticle*'s parallel to its list of straw hits
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::associateMcTruth() {

    const ComboHit* hit0 = &_data->chcol->at(0);


    for (int is=0; is<kNStations; is++) {
      int nseeds = _data->seedHolder[is].size();
      for (int i=0; i<nseeds; ++i) {
	DeltaSeed* seed = _data->seedHolder[is].at(i);
//-----------------------------------------------------------------------------
// define MC pointers (SimParticle's) for the first two, "pre-seed", hits
//-----------------------------------------------------------------------------
	if (seed->fType != 0) {
	  int loc1 = seed->fHitData[0]->fHit-hit0;
	  int loc2 = seed->fHitData[1]->fHit-hit0;
	
	  seed->fPreSeedMcPart[0] = _list_of_mc_part_hit.at(loc1);
	  seed->fPreSeedMcPart[1] = _list_of_mc_part_hit.at(loc2);
	}
	else {
	  seed->fPreSeedMcPart[0] = NULL;
	  seed->fPreSeedMcPart[1] = NULL;
	}
//-----------------------------------------------------------------------------
// define MC pointers (McPart_t's) for hits in all faces
//-----------------------------------------------------------------------------
	for (int face=0; face<kNFaces; face++) {
	  int nh = seed->hitlist[face].size();
	  for (int ih=0; ih<nh; ih++) {
	    const ComboHit* hit = seed->hitlist[face][ih]->fHit;
	    int loc = hit-hit0;
	    // get the StrawDigi indices associated with this ComboHit
	    std::vector<StrawDigiIndex> shids;
	    _data->chcol->fillStrawDigiIndices(*(_data->event),loc,shids);
	    McPart_t* mc = _list_of_mc_part_hit.at(shids[0]);//loc);
	    seed->fMcPart[face].push_back(mc);
//-----------------------------------------------------------------------------
// count CE hits 
//-----------------------------------------------------------------------------
	    if ((mc->fPdgID == 11) && (mc->fStartMom > 95) && (mc->fStartMom <110)) {
	      seed->fNHitsCE += 1;
	      if (seed->fDeltaIndex >= 0) {
		DeltaCandidate* dc = &_data->deltaCandidateHolder[seed->fDeltaIndex];
		dc->fNHitsCE += 1;
	      }
	    }
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
    int          npart;

    int ndelta = _data->deltaCandidateHolder.size();
    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = & _data->deltaCandidateHolder[idelta];
      npart              = 0;

      for (int is=dc->fFirstStation; is<=dc->fLastStation; is++) {
	DeltaSeed* ds = dc->Seed(is);
	if (ds == 0) continue;
	for (int face=0; face<kNFaces; face++) {
	  int nh = ds->hitlist[face].size();
	  for (int ih=0; ih<nh; ih++) {
//-----------------------------------------------------------------------------
// assign delta index to each hit flagged as delta
//-----------------------------------------------------------------------------
	    HitData_t* hd = (HitData_t*) ds->hitlist[face][ih];
	    hd->fDeltaIndex = idelta;
//-----------------------------------------------------------------------------
// try to identify a reconstructed delta-electron with the MC particle
//-----------------------------------------------------------------------------
	    McPart_t* mc = ds->fMcPart[face][ih];  // list parallel to the list of hits
	    
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
      }
//-----------------------------------------------------------------------------
// look at the particles contributed charge to the seed and determine 
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
    int nmc    = _list_of_mc_particles.size();

    _nDeltaHitsTot  = 0;
    _nDeltaHitsReco = 0;

    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
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

      if (mc->fPdgID == 11) {
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
	//	_listOfMcStrawHits = _mcUtils->getListOfMcStrawHits(_data->event, _stepPointMcCollTag);
	InitMcDiag();
	associateMcTruth();
      }
    }
    
    if (_printDeltaSeeds != 0) {
      for (int st=0; st<kNStations; ++st) {
	int nseeds = _data->seedHolder[st].size();
	printf("station: %2i N(seeds): %3i\n",st,nseeds);
	if (nseeds > 0) {
	  printf("---------------------------------------------------------------------------------------------------------------------------------------\n");
	  printf("      st  i  good:type   SHID:MCID(0)    SHID:MCID(1)    chi2all/N  chi21    chi22 mintime  maxtime      X        Y         Z  nfwh nht\n");
	  printf("---------------------------------------------------------------------------------------------------------------------------------------\n");
	  for (int ps=0; ps<nseeds; ++ps) {
	    DeltaSeed* seed = _data->seedHolder[st].at(ps);

	    printf("seed %2i:%03i %5i %2i ",st,ps,seed->fGood,seed->fType);
	    if (seed->fType != 0) {
	      printf("(%5i:%9i)",seed->fHitData[0]->fHit->strawId().straw()/*strawIndex().asInt()*/,seed->fPreSeedMcPart[0]->fID);//FIXME!
	      printf("(%5i:%9i)",seed->fHitData[1]->fHit->strawId().straw()/*strawIndex().asInt()*/,seed->fPreSeedMcPart[1]->fID);//FIXME!
	    }
	    else {
	      printf("(%5i:%9i)",-1,-1);
	      printf("(%5i:%9i)",-1,-1);
	    }
	    printf(" %8.2f %8.2f %8.2f",seed->Chi2AllDof(),seed->fChi21,seed->fChi22);
	    printf("%8.1f %8.1f",seed->fMinTime,seed->fMaxTime);
	    printf(" %8.3f %8.3f %9.3f",seed->CofM.x(),seed->CofM.y(),seed->CofM.z());
	    printf("%4i",seed->fNFacesWithHits);
	    printf("%4i",seed->fNHitsTot);
	    printf("\n");
//-----------------------------------------------------------------------------
// print hit ID's in each face
//-----------------------------------------------------------------------------
	    for (int face=0; face<kNFaces; face++) {
	      int first_line=1;
	      printf("       ");
	      int nh = seed->NHits(face);
	      printf("         %i:%2i ",face,nh);
	      int nprinted = 0;
	      for (int ih=0; ih<nh; ih++) {
		const ComboHit* hit =  seed->hitlist[face][ih]->fHit;

		McPart_t* mcp = seed->fMcPart[face][ih];
		int mcid = mcp->fID;
		if ((nprinted == 0) && (first_line == 0)) printf("%21s","");
		printf("(%5i:%9i)",hit->strawId().straw()/*strawIndex().asInt()*/,mcid);
		nprinted++;
		if (nprinted == 5) {
		  printf("\n");
		  first_line = 0;
		  nprinted = 0;
		}
	      }
	      if (nprinted > 0) {
		printf("\n");
	      }
	      if (nh == 0) {
		printf("(   -1:       -1)");
		printf("\n");
	      }
	    }
	  }
	}
      }
    }

//-----------------------------------------------------------------------------
// print reconstructed delta candidates
//-----------------------------------------------------------------------------
    if (_printDeltaCandidates != 0) {
      int nd = _data->deltaCandidateHolder.size();
      printf(" [DeltaFinder::debug] N(delta candidates) = %5i\n",nd);

      if (nd > 0) {
	for (int i=0; i<nd; i++) {
	  DeltaCandidate* dc = &_data->deltaCandidateHolder.at(i);
	  int pdg_id = -1;
	  if (dc->fMcPart) pdg_id = dc->fMcPart->fPdgID;
	  printf("--------------------------------------------------------------------------------------------------------------------------\n");
	  printf("      i  nh n(CE) ns s1  s2     X        Y        Z     chi21   chi22   htmin   htmax   t0min   t0max     PdgID N(MC hits)\n");
	  printf("--------------------------------------------------------------------------------------------------------------------------\n");
	  printf(":dc:%03i %3i  %3i",i,dc->fNHits,dc->fNHitsCE);
	  printf(" %3i",dc->n_seeds);
	  printf(" %2i  %2i %7.2f %7.2f %9.2f",dc->fFirstStation,dc->fLastStation,
		 dc->CofM.x(),dc->CofM.y(),dc->CofM.z());
	  printf("                            %30i %5i",pdg_id,dc->fNHitsMcP);
	  printf("\n");

	  for (int is=dc->fFirstStation;is<=dc->fLastStation; is++) {
	    DeltaSeed* ds = dc->seed[is];
	    if (ds != NULL) {
	      printf("        %3i  %3i    %3i:%03i",ds->fNHitsTot,ds->fNHitsCE,is,ds->fNumber);
	      printf(" %7.2f %7.2f %9.2f",
		     ds->CofM.x(),ds->CofM.y(),ds->CofM.z());
	      printf(" %7.1f %7.1f",ds->fChi21, ds->fChi22);
	      printf(" %7.1f %7.1f",ds->fMinTime,ds->fMaxTime);
	      printf(" %7.1f %7.1f",dc->fT0Min[is],dc->fT0Max[is]);
	      if (ds->fType != 0) {
		int f0 = ds->fType / 10;
		int f1 = ds->fType % 10;

		printf(" (%5i:%9i, %5i:%9i)",
		       ds->hitlist[f0][0]->fHit->strawId().straw()/*strawIndex().asInt()*/,ds->fMcPart[f0][0]->fID,
		       ds->hitlist[f1][0]->fHit->strawId().straw()/*strawIndex().asInt()*/,ds->fMcPart[f1][0]->fID
		       );//FIXME!
	      }
	      else {
		printf(" (%5i:%9i, %5i:%9i)",-1,-1,-1,-1);
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
      int nmc = _list_of_mc_particles.size();

      for (int i=0; i<nmc; i++) {
	McPart_t* mc = _list_of_mc_particles.at(i);

	if ((mc->fPdgID     == 11                     ) && 
	    (mc->Time()     > 550                     ) && 
	    (mc->Momentum() >  _printElectronsMinMom  ) && 
	    (mc->Momentum() <  _printElectronsMaxMom  ) && 
	    (mc->NHits()    >= _printElectronsMinNHits)    ) {

	  float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

	  if (fr < _printElectronsMaxFReco) {

	    int delta_id(-1), nseeds(0);
	    if (mc->fDelta) {
	      delta_id = mc->fDelta->fNumber;
	      nseeds   = mc->fDelta->n_seeds;
	    }

	    printf(" event: %4i electron.sim.id: %10i",_data->event->event(),mc->fID);
	    printf(" mom = %7.3f time: %8.3f deltaID: %3i nseeds: %2i nhits: %3i/%3i stations:%2i:%2i",
		   mc->Momentum(), mc->Time(), 
		   delta_id, nseeds,
		   mc->fNHitsDelta, 
		   mc->NHits(), 
		   mc->fFirstStation, mc->fLastStation);
	    printf(" freco: %6.3f\n",fr);

	    if (_printElectronsHits > 0) {
	      int nh = mc->fListOfHits.size();
	      if (nh > 0) printHitData(NULL,-1);
	      for (int ih=0; ih<nh; ih++) {
		printHitData(mc->fListOfHits[ih],ih);
	      }
	    }
	  }
	}
      }
    }

    if (_printOTracker > 0) printOTracker();

    if (_printShcol) printComboHitCollection();

    return 0;
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printHitData(const HitData_t* Hd, int Index) {

    if (Index < 0) {
      printf("#-----------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------\n");
      printf("#      SHID  St:Pl P L Str     Time     dt        eDep       wdist     wres   ");
      printf("     PDG           ID       p      X        Y         Z   DeltaID radOK edepOK\n");
      printf("#-----------------------------------------------------------------------------");
      printf("------------------------------------------------------------------------------\n");
      return;
    }

    const ComboHit* sh0 = &_data->chcol->at(0);
    const ComboHit* sh  = Hd->fHit;
    int loc             = sh-sh0;

    // const StrawHitPosition* shp = Hd->fPos;
    const StrawHitFlag*     shf = &_data->shfcol->at(loc);
   
    int radselOK        = (! shf->hasAnyProperty(StrawHitFlag::radsel));
    int edepOK          = (! shf->hasAnyProperty(StrawHitFlag::energysel));

    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999);
    float              mc_mom(-9999.);
	
    if (_mcDiag) {
      sim    = _mcUtils->getSimParticle(_data->event,loc);
      pdg_id = _mcUtils->getPdgID(sim);
      sim_id = _mcUtils->getID(sim);
      mc_mom = _mcUtils->getStartMom(sim);
    }

    // const mu2e::Straw* straw = &_data->tracker->getStraw(sh->strawIndex());
    
    printf("%5i ",loc);
    printf("%5i" ,sh->strawId().straw()/*strawIndex().asInt()*/);
	
    printf("  %2i:%2i %1i %1i %2i   %8.3f %7.3f  %9.6f   %8.3f %8.3f %10i   %10i %8.3f %8.3f %8.3f %9.3f %5i %5i %5i\n",
	   sh->strawId().station(),//straw->id().getStation(),
	   sh->strawId().plane(),  //straw->id().getPlane(),
	   sh->strawId().panel(),  //straw->id().getPanel(),
	   sh->strawId().layer(),  //straw->id().getLayer(),
	   sh->strawId().straw(),  //straw->id().getStraw(),
	   sh->time(),
	   -1.,//sh->dt(),//FIXME!
	   sh->energyDep(),
	   sh->wireDist(),
	   sh->posRes(ComboHit::wire),
	   pdg_id,
	   sim_id,
	   mc_mom,
	   sh->pos().x(),
	   sh->pos().y(),
	   sh->pos().z(),
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
	for (int ip=0; ip<kNPanelsPerFace; ip++) {
	  PanelZ_t* pz = &_data->oTracker[is][face][ip];
	  printf("#        --------------- station: %2i face: %2i panel: %2i nhits:%3li\n",
		 is,face,ip, pz->fHitData.size());
	  int nh2 = pz->fHitData.size();// pz->fHitData[0].size()+pz->fHitData[1].size();
	  if (nh2 > 0) printHitData(NULL,-1);
	  // for (int il=0; il<2; il++) {
	  int nh = pz->fHitData.size();
	  for (int ih=0; ih<nh; ih++) {
	    printHitData(&pz->fHitData[ih],ih);
	  }
	  nhitso += nh;
	  // }
	}
      }
    }

    printf(" nhits, nhitso : %6i %6i \n", (int) _data->chcol->size(),nhitso);
  }

//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printComboHit(const ComboHit* Sh, int Index) {

    if (Index <= 0) {
      printf("#-------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");
      printf("#S:F  I   SHID  Plane   Panel  Layer   Straw     Time          dt       eDep ");
      printf("           PDG         ID         p   radselOK edepOK\n");
      printf("#-------------------------------------------------------------------------");
      printf("-----------------------------------------------------\n");
      if (Index < 0) return;
    }

    const ComboHit* sh0 = &_data->chcol->at(0);
    int loc             = Sh-sh0;
    
    const StrawHitFlag* shf = &_data->shfcol->at(loc);
   
    int radselOK        = (! shf->hasAnyProperty(StrawHitFlag::radsel));
    int edepOK          = (! shf->hasAnyProperty(StrawHitFlag::energysel));

    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999);
    float              mc_mom(-9999.);
	
    if (_mcDiag) {
      sim    = _mcUtils->getSimParticle(_data->event,loc);
      pdg_id = _mcUtils->getPdgID(sim);
      sim_id = _mcUtils->getID(sim);
      mc_mom = _mcUtils->getStartMom(sim);
    }

    // const mu2e::Straw* straw = &_data->tracker->getStraw(Sh->strawIndex());
    
    printf("%5i ",loc);
    printf("%5i" ,Sh->strawId().straw());//FIXME! Sh->strawIndex().asInt());
	
    printf("  %5i  %5i   %5i   %5i   %8.3f   %8.3f   %9.6f   %10i   %10i  %8.3f %5i %5i\n",
	   Sh->strawId().plane(),  //straw->id().getPlane(),
	   Sh->strawId().panel(),  //straw->id().getPanel(),
	   Sh->strawId().layer(),  //straw->id().getLayer(),
	   Sh->strawId().straw(),  //straw->id().getStraw(),
	   Sh->time(),
	   -1.,                //Sh->dt(),//FIXME!
	   Sh->energyDep(),
	   pdg_id,
	   sim_id,
	   mc_mom,
	   radselOK,edepOK);
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printComboHitCollection() {

    int nsh = _data->chcol->size();
    printf(" ComboHitCollection: nhits : %6i\n", nsh);

    for (int i=0; i<nsh; i++) {
      const ComboHit* sh = &_data->chcol->at(i);
      printComboHit(sh,i);
    }

  }

}

DEFINE_ART_CLASS_TOOL(mu2e::DeltaFinderDiag)

#endif
