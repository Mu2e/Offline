#ifndef __CalPatRec_DeltaFinderDiag_hh__
#define __CalPatRec_DeltaFinderDiag_hh__

#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "CalPatRec/inc/DeltaFinder_types.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "CalPatRec/inc/ModuleHistToolBase.hh"
#include "CalPatRec/inc/McUtilsToolBase.hh"

using namespace std;

namespace mu2e {
  using namespace DeltaFinderTypes;

  class SimParticle;
  class PtrStepPointMCVectorCollection;
  
  class DeltaFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets =  10,
      kNSeedHistSets  =  10,
      kNDeltaHistSets =  10,
      kNMcHistSets    = 200
    };

    struct SeedHist_t {
      TH1F*  fChi2Dof;
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
    int                                   _printDeltaSeeds;
    int                                   _printDeltaCandidates;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    int                                   _firstCall;

    const PtrStepPointMCVectorCollection* _listOfMcStrawHits;
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
    
    void        printStrawHit(const StrawHit* Sh, int Index);
    void        printOTracker();
//-----------------------------------------------------------------------------
// overriden virtual functions of the base class
//-----------------------------------------------------------------------------
  public:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1 ) override ;
    virtual int debug         (void* Data) override ;
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
    _printDeltaSeeds       (PSet.get<int>          ("printDeltaSeeds"              )),
    _printDeltaCandidates  (PSet.get<int>          ("printDeltaCandidates"         ))
  {
    printf(" DeltaFinderDiag::DeltaFinderDiag : HOORAY! \n");
    _firstCall   =  1;
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

    Hist->fChi2Dof         = Dir->make<TH1F>("chi2dof"     , "Chi squared/degrees of freedom", 1000, 0., 100.);
    Hist->fNFacesWithHits  = Dir->make<TH1F>("nfaces_wh"   , "Number of faces with hits", 5, 0., 5.);
    Hist->fNHitsPerFace    = Dir->make<TH1F>("nhits_face"  , "Number of hits per face", 20, 0., 20.);
    Hist->fNHitsPerSeed    = Dir->make<TH1F>("nhits_seed"  , "Number of hits per seed", 40, 0., 40.);
    Hist->fSeedRadius      = Dir->make<TH1F>("rad"         , "Seed radius", 500, 0., 1000.);
    Hist->fSeedMomentum    = Dir->make<TH1F>("mom"         , "Seed momentum", 400, 0., 400.);
    Hist->fSeedSize        = Dir->make<TH2F>("preseed_size", "Seed (nh2+1):(nh1+1)", 20, 0., 20.,20,0,20);
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
    Hist->fMaxSeg     = Dir->make<TH1F>("max_seg", "Max N Segments", 20, 0., 20.);

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

	  for (int l=0; l<2; ++l) { //loop over layers
	    int nhits =  panelz->fHitData[l].size();
	    for (int i=0; i<nhits; i++) {
	      int counter  = panelz->fHitData[l].at(i).fNSecondHits;
	      Hist->fNSecondHits->Fill(counter);
	    }
	  }
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

    Hist->fChi2Dof->Fill(Seed->chi2dof);

    int loc1 = Seed->fType/10;
    int loc2 = Seed->fType % 10;

    int nh1 = Seed->hitlist[loc1].size();
    int nh2 = Seed->hitlist[loc2].size();

    Hist->fSeedSize->Fill(nh1+1,nh2+1);
    Hist->fNFacesWithHits->Fill(Seed->fNFacesWithHits);
    Hist->fSeedRadius->Fill(Seed->CofM.perp());

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
    if (_mcDiag) {
      _listOfMcStrawHits = _mcUtils->getListOfMcStrawHits(_data->event, _stepPointMcCollTag);

      InitMcDiag();
      associateMcTruth();
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
	DeltaSeed* seed = &_data->seedHolder[s].at(se);
	
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
      //      const SimParticle* sim = mc->fSim;

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

	    if (_data->debugLevel > 0) {
//-----------------------------------------------------------------------------
// a closer look at misreconstructed delta electrons
//-----------------------------------------------------------------------------
	      float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

	      if ((mc->Momentum() < 5) && (mc->Time() > 550) && (mc->NHits() > 40) && (fr < 0.5)) {
		printf(" event: %6i missed delta: sim.id = %10i mom = %10.3f time= %9.3f nhits = %3i nhits(delta): %3i first: %2i last: %2i",
		       _data->event->event(),
		       mc->fID, mc->Momentum(), mc->Time(), 
		       mc->NHits(), mc->fNHitsDelta, 
		       mc->fFirstStation, mc->fLastStation);
		printf(" fraction: %6.3f\n",fr);
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

    int nsh = _data->shcol->size();
    _list_of_mc_part_hit.resize(nsh);

    for (int ish=0; ish<nsh; ish++) {
      
      const StrawHit* sh = &_data->shcol->at(ish);
      const mu2e::SimParticle* sim = _mcUtils->getSimParticle(_listOfMcStrawHits,ish);
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
      
      mc->fListOfHits.push_back(sh);

      StrawIndex shid = sh->strawIndex();
      const Straw& straw  = _data->tracker->getStraw(shid);
      int station = straw.id().getStation();
      if (station < mc->fFirstStation) mc->fFirstStation = station;
      if (station > mc->fLastStation ) mc->fLastStation  = station;

      if (sh->time() < mc->fTime) mc->fTime = sh->time();

      _list_of_mc_part_hit[ish] = mc;
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

    const StrawHit* hit0 = &_data->shcol->at(0);

    for (int is=0; is<kNStations; is++) {
      int nseeds = _data->seedHolder[is].size();
      for (int i=0; i<nseeds; ++i) {
	DeltaSeed* seed = &_data->seedHolder[is].at(i);
//-----------------------------------------------------------------------------
// define MC pointers (SimParticle's) for the first two, "pre-seed", hits
//-----------------------------------------------------------------------------
	int loc1 = seed->fHit[0]-hit0;
	int loc2 = seed->fHit[1]-hit0;
	
	seed->fPreSeedMcPart[0] = _list_of_mc_part_hit.at(loc1);
	seed->fPreSeedMcPart[1] = _list_of_mc_part_hit.at(loc2);
//-----------------------------------------------------------------------------
// define MC pointers (McPart_t's) for hits in all faces
//-----------------------------------------------------------------------------
	for (int face=0; face<kNFaces; face++) {
	  int nh = seed->hitlist[face].size();
	  for (int ih=0; ih<nh; ih++) {
	    const StrawHit* hit = seed->hitlist[face].at(ih);
	    int loc = hit-hit0;
	    seed->fMcPart[face].push_back(_list_of_mc_part_hit.at(loc));
	  }
	}
      }
    }
//-----------------------------------------------------------------------------
// for each MC electron calculate the number of reconstructed hits
//-----------------------------------------------------------------------------
    int nmc    = _list_of_mc_particles.size();
    int ndelta = _data->deltaCandidateHolder.size();

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
	const StrawHit* hit = mc->fListOfHits.at(ih);
	int hit_is_found = 0;
	for (int id=0; id<ndelta; id++) {
	  DeltaCandidate* dc = &_data->deltaCandidateHolder.at(id);
//-----------------------------------------------------------------------------
// now loop over the DeltaCandidate hits - the idea is that all of them will be marked 
// as 'delta electron' hits and not used in the pattern recognition
//-----------------------------------------------------------------------------
	  for (int is=dc->st_start; is<=dc->st_end; is++) {
	    DeltaSeed* ds = dc->seed[is];
	    if (ds) {
	      for (int iface=0; iface<kNFaces; iface++) {
		int nh2 = ds->NHits(iface);
		for (int ih2=0; ih2<nh2; ih2++) {
		  const StrawHit* delta_hit = ds->hitlist[iface][ih2];
	
		  if (delta_hit == hit) { 
		    hit_is_found = 1;
		    break;
		  }
		}
	      }
	    }
	  }
	  if (hit_is_found) break;
	}

	if (hit_is_found) mc->fNHitsDelta += 1;
      }

      if (mc->fPdgID == 11) {
	float mom = mc->Momentum();
	if (mom < 20) {
	  _nDeltaHitsTot  += mc->NHits();
	  _nDeltaHitsReco += mc->fNHitsDelta;
	}
      }
    }
//-----------------------------------------------------------------------------
// proceed with the delta candidates
//-----------------------------------------------------------------------------
    int const    max_part(1000);
    McPart_t*    part [max_part];
    int          nhits[max_part];
    int          npart;

    for (int i=0; i<ndelta; i++) {
      DeltaCandidate* d = & _data->deltaCandidateHolder.at(i);
      npart             = 0;

      for (int is=d->st_start; is<=d->st_end; is++) {
	if (! d->st_used[is]) continue;
	DeltaSeed* s = d->seed[is];
	for (int face=0; face<kNFaces; face++) {
	  int nh = s->hitlist[face].size();
	  for (int ih=0; ih<nh; ih++) {
	    McPart_t* sim = s->fMcPart[face][ih];

	    int found = 0;

	    for (int ip=0; ip<npart; ip++) {
	      if (sim == part[ip]) {
		found = 1;
		nhits[ip]++;
		break;
	      }
	    }

	    if (found == 0) {
	      if (npart < max_part) {
		part [npart] = sim;
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
	d->fMcPart         = part [ibest];
	d->fNHitsMcP       = nhits[ibest];
	d->fMcPart->fDelta = d;
      }
    }
    return 0;
  }
  
//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  int DeltaFinderDiag::debug(void* Data) {
//-----------------------------------------------------------------------------
// print DeltaSeeds - pieces of delta electrons reconstructed within one station
//-----------------------------------------------------------------------------
    _data = (Data_t*) Data;
    
    if (_data->printDeltaSeeds != 0) {
      for (int st=0; st<kNStations; ++st) {
	int nseeds = _data->seedHolder[st].size();
	printf("station: %2i N(seeds): %3i\n",st,nseeds);
	if (nseeds > 0) {
	  printf("------------------------------------------------------------------------------------------------------\n");
	  printf("      st  i  good type   SHID:MCID(0)    SHID:MCID(1)       chi2  mintime maxtime   X        Y         Z    nfwh  nht\n");
	  printf("------------------------------------------------------------------------------------------------------\n");
	  for (int ps=0; ps<nseeds; ++ps) {
	    DeltaSeed* seed = &_data->seedHolder[st].at(ps);

	    printf("seed %3i %2i %4i %2i ",st,ps,seed->fGood,seed->fType);
	    printf("(%5i:%9i)",seed->fHit[0]->strawIndex().asInt(),seed->fPreSeedMcPart[0]->fID);
	    printf("(%5i:%9i)",seed->fHit[1]->strawIndex().asInt(),seed->fPreSeedMcPart[1]->fID);
	    printf(" %8.2f",seed->chi2dof);
	    printf("%8.1f %8.1f",seed->fMinTime,seed->fMaxTime);
	    printf(" %8.3f %8.3f %9.3f",seed->CofM.x(),seed->CofM.y(),seed->CofM.z());
	    printf("%4i",seed->fNFacesWithHits);
	    printf("%4i",seed->fNHitsTot);
	    printf("\n");
	    //-----------------------------------------------------------------------------
	    // print hit ID's in each face
	    //-----------------------------------------------------------------------------
	    for (int face=0; face<kNFaces; face++) {
	      printf("       ");
	      int nh = seed->NHits(face);
	      printf("       %i:%3i ",face,nh);
	      for (int ih=0; ih<nh; ih++) {
		const StrawHit* hit =  seed->hitlist[face][ih];

		McPart_t* mcp = seed->fMcPart[face][ih];
		int mcid = mcp->fID;
		printf("(%5i:%9i)",hit->strawIndex().asInt(),mcid);
	      }
	      if (nh == 0) printf("(   -1:       -1)");
	      printf("\n");
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
	  printf("----------------------------------------------------------------------\n");
	  printf("   i  nh  ns s1 s2       X         Y         Z    mintime  maxtime    \n");
	  printf("----------------------------------------------------------------------\n");
	  printf("%3i: %3i",i,dc->fNHits);
	  printf(" %3i",dc->n_seeds );
	  printf(" %2i",dc->st_start);
	  printf(" %2i",dc->st_end  );
	  printf("  %9.2f %9.2f %9.2f",dc->CofM.x(),dc->CofM.y(),dc->CofM.z());
	  printf("\n");

	  for (int is=dc->st_start;is<=dc->st_end; is++) {
	    DeltaSeed* ds = dc->seed[is];
	    if (ds != NULL) {
	      printf("     %3i     %2i:%3i",ds->fNHitsTot,is,ds->fNumber);
	      printf(" %9.2f %9.2f %9.2f",
		     ds->CofM.x(),ds->CofM.y(),ds->CofM.z());
	      printf( "%8.1f %8.1f",ds->fMinTime,ds->fMaxTime);
	      printf(" (%6i %6i)",ds->fHit[0]->strawIndex().asInt(),ds->fHit[1]->strawIndex().asInt());

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

	if ((mc->fPdgID    == 11                   ) && (mc->Time()   > 550                     ) && 
	    (mc->Momentum() > _printElectronsMinMom) && (mc->NHits()  >= _printElectronsMinNHits)    ) {

	  float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

	  if (fr < _printElectronsMaxFReco) {

	    int delta_id(-1), nseeds(0);
	    if (mc->fDelta) {
	      delta_id = mc->fDelta->fNumber;
	      nseeds   = mc->fDelta->n_seeds;
	    }

	    printf(" event: %3i electron: sim.id = %10i",_data->event->event(),mc->fID);
	    printf(" mom = %7.3f time: %8.3f nhits: %3i deltaID: %3i nseeds: %2i nhits(delta): %3i stations:%2i:%2i",
		   mc->Momentum(), mc->Time(), 
		   mc->NHits(), 
		   delta_id, nseeds,
		   mc->fNHitsDelta, 
		   mc->fFirstStation, mc->fLastStation);
	    printf(" freco: %6.3f\n",fr);

	    if (_printElectronsHits > 1) {
	      int nh = mc->fListOfHits.size();
	      for (int ih=0; ih<nh; ih++) {
		printStrawHit(mc->fListOfHits.at(ih),ih);
	      }
	    }
	  }
	}
      }
    }

    if (_printOTracker > 0) printOTracker();

    return 0;
  }
    
//-----------------------------------------------------------------------------
  void DeltaFinderDiag::printStrawHit(const StrawHit* Sh, int Index) {

    if (Index <= 0) {
      printf("--------------------------------------------------------------------------");
      printf("--------------------------------------\n");
      printf("   I   SHID  Plane   Panel  Layer   Straw     Time          dt       eDep ");
      printf("           PDG         ID         p   \n");
      printf("--------------------------------------------------------------------------");
      printf("--------------------------------------\n");
      return;
    }

    const StrawHit* sh0 = &_data->shcol->at(0);
    int loc             = Sh-sh0;
    
    const SimParticle* sim(0);
    int                pdg_id(-9999), sim_id(-9999);
    float              mc_mom(-9999.);
	
    if (_mcDiag) {
      sim    = _mcUtils->getSimParticle(_listOfMcStrawHits,loc);
      pdg_id = _mcUtils->getPdgID(sim);
      sim_id = _mcUtils->getID(sim);
      mc_mom = _mcUtils->getStartMom(sim);
    }

    const mu2e::Straw* straw = &_data->tracker->getStraw(Sh->strawIndex());
    
    printf("%5i ",loc);
    printf("%5i" ,Sh->strawIndex().asInt());
	
    printf("  %5i  %5i   %5i   %5i   %8.3f   %8.3f   %9.6f   %10i   %10i  %8.3f\n",
	   straw->id().getPlane(),
	   straw->id().getPanel(),
	   straw->id().getLayer(),
	   straw->id().getStraw(),
	   Sh->time(),
	   Sh->dt(),
	   Sh->energyDep(),
	   pdg_id,
	   sim_id,
	   mc_mom);
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
	  printf("        --------------- station: %2i face: %2i panel: %2i nhits[0]:%3li nhits[1]:%3li \n",
		 is,face,ip, pz->fHitData[0].size(),pz->fHitData[1].size());

	  for (int il=0; il<2; il++) {
	    int nh = pz->fHitData[il].size();
	    for (int ih=0; ih<nh; ih++) {
	      printStrawHit(pz->fHitData[il].at(ih).fHit,ih);
	    }
	    nhitso += nh;
	  }
	}
      }
    }

    printf(" nhits, nhitso : %6i %6i \n", (int) _data->shcol->size(),nhitso);
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::DeltaFinderDiag)

#endif
