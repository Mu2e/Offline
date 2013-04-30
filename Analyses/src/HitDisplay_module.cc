//
// A sandbox for playing with tracks, including transformations to different representations.
// This is not production code but feel free to look at it.
//
// $Id: HitDisplay_module.cc,v 1.19 2013/04/30 20:29:58 murat Exp $
// $Author: murat $
// $Date: 2013/04/30 20:29:58 $
//
// Original author Rob Kutschke.
//
// small black triangles: MC truth on the trajectory
// red hits: hits on a track reconstructed as a downstream-moving conversion electron 
// 

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"

#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkHotList.hh"
#include "KalmanTrack/KalHit.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"

// ROOT includes
#include "TApplication.h"
#include "TArc.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TBox.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TText.h"
#include "TNtuple.h"

// Other includes
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class HitDisplay : public art::EDAnalyzer {
  private:

    // The module label of this instance of this module.
    std::string moduleLabel_;

    // Label of the modules that created the data products.
    std::string generatorModuleLabel_;
    std::string g4ModuleLabel_;
    std::string hitMakerModuleLabel_;

    // Name of the tracker StepPoint collection
    std::string trackerStepPoints_;

    // Cuts used inside SimParticleWithHits:
    //  - drop hits with too little energy deposited.
    //  - drop SimParticles with too few hits.
    double minEnergyDep_;
    double timeWindow_;
    size_t minHits_;

    // Options to control the display
    bool doDisplay_;
    bool clickToAdvance_;
    bool printHits_;

    TApplication* application_;
    TDirectory*   directory_;
    TCanvas*      canvas_;

    TH1F* _hnHits;
    TH1F* _hEnergyDep;
    TH1F* _hDeltaT;
    TH1F* _hTime;
    TH1F* _hx;
    TH1F* _hxnorm;

    TH1F* _hMissDist;

    //     TNtuple* _ntTrack;
    //     TNtuple* _ntHit;

    const VaneCalorimeter* fCal;

    int       fNClusters;
    // 4 hypotheses: dem, uep, dmm, ump
    int       fNTracks[4];

    TMarker*  fMarker;

  public:
    explicit HitDisplay(fhicl::ParameterSet const& pset);
    virtual ~HitDisplay();

    void beginJob( );
    void analyze( art::Event const& e );
    void printCaloCluster(const CaloCluster* Cl, const char* Opt) ;
    void printKalRep(const KalRep* Trk, const char* Opt); 

  };

  HitDisplay::~HitDisplay() { 
    delete application_; 
  }

  HitDisplay::HitDisplay(fhicl::ParameterSet const& pset):
    moduleLabel_(pset.get<string>("module_label")),
    generatorModuleLabel_(pset.get<std::string>("generatorModuleLabel")),
    g4ModuleLabel_(pset.get<std::string>("g4ModuleLabel")),
    hitMakerModuleLabel_(pset.get<std::string>("hitMakerModuleLabel")),
    trackerStepPoints_(pset.get<std::string>("trackerStepPoints")),
    minEnergyDep_(pset.get<double>("minEnergyDep")),
    timeWindow_(pset.get<double>("timeWindow")),
    minHits_(pset.get<unsigned>("minHits")),
    doDisplay_(pset.get<bool>("doDisplay",true)),
    clickToAdvance_(pset.get<bool>("clickToAdvance",false)),
    printHits_(pset.get<bool>("printHits",false)),
    application_(0),
    directory_(0),
    canvas_(0),
    _hnHits(0),
    _hEnergyDep(0),
    _hDeltaT(0),
    _hTime(0),
    _hx(0),
    _hxnorm(0),
    _hMissDist(0)
  {
    fCal = 0;
  }


  //-----------------------------------------------------------------------------
  void HitDisplay::beginJob(){

    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    // Create some 1D histograms.
    _hnHits     =  tfs->make<TH1F>( "hnHits",      "StrawHits per Event",              100,    0.,    200. );
    _hEnergyDep =  tfs->make<TH1F>( "hEnergyDep",  "Energy Deposition per Hit;(keV)",  100,    0.,     20. );
    _hDeltaT    =  tfs->make<TH1F>( "hDeltaT",     "Delta(time);(ns)",                 100,   -5.,      5. );
    _hTime      =  tfs->make<TH1F>( "hTime",       "Time;(ns)",                        100,    0.,   2000. );
    _hx         =  tfs->make<TH1F>( "hx",          "Displacement;(mm)",                100, -1000.,  1000. );
    _hxnorm     =  tfs->make<TH1F>( "hxnorm",      "Displacement;(HalfLength)",        100,    -1.,     1. );
    _hMissDist  =  tfs->make<TH1F>( "hMissDist",   "Distance to wire",                 100,   -50.,    50. );
    //     _ntTrack    =  tfs->make<TNtuple>( "ntTrack",  "TrackInfo", "d0gen:d0first" );
    //     _ntHit      =  tfs->make<TNtuple>( "ntHit",    "HitInfo",   "dca:z" );

    if ( !doDisplay_ ) return;

    // If needed, create the ROOT interactive environment. See note 1.
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
    }

    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + moduleLabel_;
    TString title = "Canvas for " + moduleLabel_;
    int window_size(800);
    canvas_ = tfs->make<TCanvas>(name,title,window_size,window_size);

    directory_ = gDirectory;

    fMarker = new TMarker(0,0,20);
    fMarker->SetMarkerSize(0.3);
  }


  //-----------------------------------------------------------------------------
  // 'banner' : print banner
  // 'data'   : print track data
  // 'hits'   : print hits
  void HitDisplay::printKalRep(const KalRep* Trk, const char* Opt) {

    TString opt = Opt;

    if ((opt == "") || (opt == "banner")) {
      printf("------------------------------------------------------------------------------------------------\n");
      printf("TrkID  Q   momentum       pt       costh      d0        z0         T0     Nact     chi2 N(dof) FitConsistency\n");
      printf("------------------------------------------------------------------------------------------------\n");
    }
    
    if ((opt == "") || (opt.Index("data") >= 0)) {
      //      Trk->printAll();
      double mom   = Trk->momentum().mag();
      double pt    = Trk->momentum().perp();
      double costh = Trk->momentum().cosTheta();
      double d0    = Trk->helix(0).d0();
      double z0    = Trk->helix(0).z0();
      double chi2  = Trk->chisq();
      int    ndof  = Trk->nDof ();
      int    nact  = Trk->nActive();
      double t0    = Trk->t0().t0();
      double fit_consistency = Trk->chisqConsistency().consistency();
      int    q     = Trk->charge();

      printf("%2i %10.3f %10.3f %8.3f %8.3f %10.3f %10.3f %3i %10.3f %3i %10.3f\n",
	     q,mom,pt,costh,d0,z0,t0,nact,chi2,ndof,fit_consistency);

      if (opt.Index("hits") >= 0) {
	//-----------------------------------------------------------------------------
	// print detailed information about the track hits
	//-----------------------------------------------------------------------------
	const TrkHotList* hot_list = Trk->hotList();

	printf("---------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
	printf(" ih U A     len      rms       x          y          z       HitT     HitDt");
	printf("  SInd  Dev Sec Lay  N  Iamb     T0    Rdrift     Xs         Ys          Zs        resid\n");
	printf("---------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");

	int i = 0;
 	for(TrkHotList::hot_iterator it=hot_list->begin(); it<hot_list->end(); it++) {
	  //	  const KalHit* kh = (*it).kalHit();

	  // TrkStrawHit inherits from TrkHitOnTrk

  	  const TrkStrawHit* hit = (const TrkStrawHit*) &(*it);

	  const StrawHit* sh = &hit->strawHit();
	  Straw*   straw = (Straw*) &hit->straw();

	  double len = hit->fltLen();

	  HepPoint  plen = Trk->position(len);

	  printf("%3i %1i %1i %10.3f %6.3f %10.3f %10.3f %10.3f %8.3f %7.3f",
		 ++i,
		 hit->isUsable(),
		 hit->isActive(),
		 len,
		 hit->hitRms(),
		 plen.x(),plen.y(),plen.z(),
		 sh->time(), sh->dt()
		 );

	  Hep3Vector pos;
	  hit->hitPosition(pos);
	  printf("%6i %3i %3i %3i %3i %3i %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f\n",
		 straw->index().asInt(), 
		 straw->id().getDevice(),
		 straw->id().getSector(),
		 straw->id().getLayer(),
		 straw->id().getStraw(),
		 
		 hit->ambig(),
		 hit->hitT0().t0(),
		 hit->driftRadius(),
		 pos.x(),
		 pos.y(),
		 pos.z(),
		 hit->resid()
		 );
	  //  	  trkhit->print(std::cout);
 	}
      }
    }
  }

  //-----------------------------------------------------------------------------
  // Opt: 
  // 'banner' : print banner
  //-----------------------------------------------------------------------------
  void HitDisplay::printCaloCluster(const CaloCluster* Cl, const char* Opt) {
    const char oname[] = "HitDisplay::printCaloCluster";

    int row, col, ir, iz;
    TString opt = Opt;

    art::ServiceHandle<GeometryService> geom;
    //if(! geom->hasElement<VaneCalorimeter>() ) return;
    GeomHandle<Calorimeter> cg;

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("-----------------------------------------------------------------------\n");
      printf("VaneID       Time   Row   Col   Energy       X          Y        Z     \n");
      printf("-----------------------------------------------------------------------\n");
    }
    
    if ((opt == "") || (opt.Index("data") >= 0)) {
      row = Cl->cogRow();
      col = Cl->cogColumn();

      if ((row < 0) || (row > 9999)) row = -9999;
      if ((col < 0) || (col > 9999)) col = -9999;

      printf("%6i  %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f\n",
	     Cl->vaneId(),
	     Cl->time(),
	     row, col,
	     Cl->energyDep(),
	     Cl->cog3Vector().x(),
	     Cl->cog3Vector().y(),
	     Cl->cog3Vector().z() );

      const CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();

      int nh = caloClusterHits.size();
      Hep3Vector pos;
      //-----------------------------------------------------------------------------
      // print individual crystals in local vane coordinate system
      //-----------------------------------------------------------------------------
      for (int i=0; i<nh; i++) {
	const CaloCrystalHit* hit = &(*caloClusterHits.at(i));
	int id = hit->id();

	pos = cg->localCrystalOrigin(id);

	if(geom->hasElement<mu2e::VaneCalorimeter>() ){
	  GeomHandle<VaneCalorimeter> cgvane;
	  iz  = cgvane->nCrystalZ();
	  ir  = cgvane->nCrystalR();
	}else{
	  HexMap hmap;
	  //	  HexLK lk = hmap.lk(id);

	  printf(">>> ERROR in %s: yet to learn how to print the crystal indices in case of disks\n", oname);

// 	  iz = lk._l;
// 	  ir = lk._k;

 	  iz = -1;
	  ir = -1;
	}

	
	
	printf("%6i  %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	       id,
	       hit->time(),
	       iz,ir,
	       hit->energyDep(),
	       pos.x(),
	       pos.y(),
	       pos.z(),
	       hit->energyDepTotal()
	       );
      }
    }
  }


  //-----------------------------------------------------------------------------
  void HitDisplay::analyze(art::Event const& event) {
    const char* name = "HitDisplay::analyze";

    // Tracker geometry.
    GeomHandle<TTracker> ttracker;

    // if (fCal == 0) {
//       GeomHandle<VaneCalorimeter> cg;
//       fCal = &(*cg);
    // }

    TText          t;
    TEllipse*      e;
    TGraph*        graph (0);

    TrackTool      *tt(0), *tg(0), *tmid(0);
    const GenParticle* gen;
 
    TObjArray      list_of_ellipses;
    int            n_displayed_hits;
    vector<double> xStep,yStep;

    CaloCluster*    cl;
    //    int             vane_id;
    double xl, yl,  event_time;
    const KalRep*   trk;
    TString         opt; 

    printf("[%-30s] RUN: %10i EVENT: %10i\n",name,event.run(),event.event());

    // Geometry of tracker envelope.
    TubsParams envelope(ttracker->getInnerTrackerEnvelopeParams());

    // Tracker calibration object.
    ConditionsHandle<TrackerCalibrations> trackCal("ignored");

    // Get information from the event.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(generatorModuleLabel_,gensHandle);
    GenParticleCollection const& gens = *gensHandle;

    //     art::Handle<SimParticleCollection> simsHandle;
    //     event.getByLabel(g4ModuleLabel_,simsHandle);
    //     SimParticleCollection const& sims = *simsHandle;

    art::Handle<StrawHitCollection> hitsHandle;
    event.getByLabel(hitMakerModuleLabel_,hitsHandle);
    StrawHitCollection const& hits = *hitsHandle;

    art::Handle<CaloClusterCollection> calo_cluster_handle;
    event.getByLabel("makeCaloCluster","AlgoCLOSESTSeededByENERGY",calo_cluster_handle);
    CaloClusterCollection* fListOfClusters;
    fListOfClusters = (CaloClusterCollection*) &(*calo_cluster_handle);
    fNClusters      = fListOfClusters->size();

    art::Handle<KalRepCollection> demHandle;
    event.getByLabel("trkPatRec1","DownstreameMinus", demHandle);
    const KalRepCollection*  dem = &(*demHandle);
    fNTracks[0] = dem->size();
    //-----------------------------------------------------------------------------
    // three other hit collections
    //-----------------------------------------------------------------------------
    art::Handle<KalRepCollection> uepHandle;
    event.getByLabel("trkPatRec2","UpstreamePlus", uepHandle);
    const KalRepCollection*  uep = &(*uepHandle);
    fNTracks[1] = uep->size();

    art::ServiceHandle<mu2e::GeometryService> geom;
    //if(! geom->hasElement<mu2e::VaneCalorimeter>() ) return;
    //mu2e::GeomHandle<mu2e::VaneCalorimeter> vc;

    //     art::Handle<KalRepCollection> dmmHandle;
    //     event.getByLabel("trkPatRec3","DownstreammuMinus", dmmHandle);
    //     const KalRepCollection*  dmm = &(*dmmHandle);
    //     fNTracks[2] = dmm->size();

    //     art::Handle<KalRepCollection> umpHandle;
    //     event.getByLabel("trkPatRec4","UpstreammuPlus", umpHandle);
    //     const KalRepCollection*  ump = &(*umpHandle);
    //     fNTracks[3] = ump->size();

    /*
      art::Handle<StrawHitMCTruthCollection> hitsTruthHandle;
      event.getByLabel(hitMakerModuleLabel_,hitsTruthHandle);
      StrawHitMCTruthCollection const& hitsTruth = *hitsTruthHandle;
    */

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(hitMakerModuleLabel_,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const& hits_mcptr = *mcptrHandle;

    art::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel(g4ModuleLabel_,trackerStepPoints_,stepsHandle);
    StepPointMCCollection const& steps = *stepsHandle;

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits simsInfo( event,
                                   g4ModuleLabel_,
                                   hitMakerModuleLabel_,
                                   trackerStepPoints_,
                                   minEnergyDep_,
                                   minHits_ );

    SimParticleCollection::key_type key(1);

    SimParticleInfo const* info(simsInfo.findOrNull(key));
    const StepPointMC* firstStep; 

    if ( !info ) {
      printf("[%-30s] *** ERROR: SimParticleInfo missing, SKIPPING EVENT: %10i\n",name,event.event());
      firstStep = 0;
      //      return;
    }
    else {
      firstStep = &info->firstStepPointMCinTracker();
    }
    
    static TLine*  line(0);
    static TArc*   arc(0);
    static TArc*   arccalo(0);
    static TArrow* arrow(0);
    static TBox*   box(0);

    if (line == 0) {
      line   = new TLine;
      arc    = new TArc;
      arccalo = new TArc;
      arrow  = new TArrow;
      box    = new TBox;

      box->SetFillStyle(3002);
      box->SetFillColor(4);
      box->SetLineColor(4);
    }

    SortedStepPoints sortedSteps(key,steps);
    
    const StepPointMC* midStep = & sortedSteps.middleByZ();
    
    if ((midStep == 0) || (firstStep ==0)) {
      printf("[%-30s] **** ERROR : HitDisplay_module::analyze : firstStep = %8p midstep = %8p, BAIL OUT\n",
	     name,firstStep, midStep);
      //      goto END_OF_ROUTINE;
    }

    Hep3Vector pos1, mom1, pos2, mom2;
    if (firstStep ) {
      pos1 = firstStep->position();
      mom1 = firstStep->momentum();
    }
    else {
      pos1.set(1.,1.,1.);
      mom1.set(1.,1.,1.);
    }

    if (midStep) {
      pos2 = midStep->position();
      mom2 = midStep->momentum();
    }
    else {
      pos2.set(1.,1.,1.);
      mom2.set(1.,1.,1.);
    }

    // The generated particle.

    gen = &gens.at(0);

    tt   = new TrackTool(gen->pdgId(), -1.,pos1,mom1,1.,Hep3Vector());
    tg   = new TrackTool(gen->pdgId(), -1.,gen->position()      ,gen->momentum()      ,1.,Hep3Vector());
    tmid = new TrackTool(gen->pdgId(), -1.,pos2  ,mom2 ,1.,Hep3Vector());

    //   _hnHits->Fill( hits.size() );


    for ( size_t ipt=0; ipt<steps.size(); ++ipt){
      StepPointMC const& step =  steps.at(ipt);
      if ( step.totalEDep() > minEnergyDep_ ) {
        xStep.push_back( step.position().x() );
        yStep.push_back( step.position().y() );
      }
    }

    arc->SetFillStyle(0);
    arccalo->SetFillStyle(0);

    //  DISPLAY:;

    if ( doDisplay_ ) {
      
      canvas_->cd(0);
      canvas_->Clear();
      //-----------------------------------------------------------------------------
      // Draw the frame
      //-----------------------------------------------------------------------------
      double plotLimits(850.);
      canvas_->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);
      
      t.SetText(-800.,900.,Form("RUN: %10i EVENT: %10i NTRACKS: %4i NCLUSTERS: %4i",
				event.run(),event.event(),fNTracks[0],fNClusters));
      t.SetTextSize(0.02);
      t.Draw();
      
      // Draw the inner and outer arc of the tracker.
      arc->SetLineColor(kBlack);
      arc->DrawArc(0.,0., envelope.outerRadius());
      arc->DrawArc(0.,0., envelope.innerRadius());
      arc->SetLineColor(kRed);
      arc->DrawArc( tt->xc(), tt->yc(), tt->rho());
      arc->SetLineColor(kMagenta);
      arc->DrawArc( tmid->xc(), tmid->yc(), tmid->rho());
      arc->SetLineColor(kRed);
      //-----------------------------------------------------------------------------
      // draw vanes or disks
      //-----------------------------------------------------------------------------
      int nv(0);// = vc->nVane();
      if( geom->hasElement<mu2e::VaneCalorimeter>() ){
	mu2e::GeomHandle<mu2e::VaneCalorimeter> vc;
	nv = vc->nVane();
      }else if( geom->hasElement<mu2e::DiskCalorimeter>() ){
	mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
	nv = dc->nDisk();
      }
      //mu2e::GeomHandle<mu2e::VaneCalorimeter> vc;
      double rmin, rmax, x1,y1,x2,y2;
      if( geom->hasElement<mu2e::VaneCalorimeter>() ){
	mu2e::GeomHandle<mu2e::VaneCalorimeter> vc;
	for (int iv=0; iv<nv; iv++) {
	  const Vane* vane = &vc->vane(iv);
	  rmin = vane->innerRadius();
	  rmax = vane->outerRadius();
	  if(iv == 3) {			// top
	    x1 = -50;
	    y1 = rmin;
	    x2 = 50;
	    y2 = rmax;
	  }
	  else if (iv == 0) { // left
	    x1 = -rmax;;
	    y1 = -50;
	    x2 = -rmin;
	    y2 = +50;
	  }
	  else if (iv == 1) {		// bottom
	    x1 = -50;
	    y1 = -rmax;
	    x2 = +50;
	    y2 = -rmin;
	  }
	  else if (iv == 2) {		// right
	    x1 = rmin;
	    y1 = -50;
	    x2 = rmax;
	    y2 = +50;
	  }
	  box->DrawBox(x1,y1,x2,y2);
	}
      }else if( geom->hasElement<mu2e::DiskCalorimeter>() ){
	
	mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
	// Draw the inner and outer arc of calorimeter

	for (int iv=0; iv<nv; iv++) {
	  
	  rmin = dc->disk(iv).innerRadius();
	  rmax = dc->disk(iv).outerRadius();
	  arccalo->SetLineColor(kGreen+iv);
	  arccalo->DrawArc(0.,0., rmax);
	  arccalo->DrawArc(0.,0., rmin);
	 
	}
      }
    }
    //-----------------------------------------------------------------------------
    // print reconstructed tracks - all 4 hypotheses for comparison...
    //-----------------------------------------------------------------------------
    printf(" [%s] NTRACKS = %4i, NCLUSTERS = %4i\n",name,fNTracks[0],fNClusters);
    printf("\n");

    opt = "data";
    if (printHits_) opt += "+hits";

    if (fNTracks[0] > 0) {
      printKalRep(0,"banner");
      for (int i=0; i<fNTracks[0]; i++ ) {
	trk = (*dem)[i];
	printf(" %2i dem ",i);
	printKalRep(trk,opt);
      }
    }

    //       for (int i=0; i<fNTracks[1]; i++ ) {
    // 	trk = (*uep)[i];
    // 	printf(" %2i uep ",i);
    // 	printKalRep(trk,opt);
    //       }

    int module_color[2] = {kRed, kMagenta};
    int color;

    if (fNClusters > 0) {
      printCaloCluster(0,"banner");

      for (int i=0; i<fNClusters; i++) {
	cl = &fListOfClusters->at(i);
	printCaloCluster(cl,opt);
	// event time defined by the most energetic cluster
	if (i == 0) {
	  event_time = cl->time()+15.;
	}
	//-----------------------------------------------------------------------------
	// display only clusters with E > 5 MeV
	//-----------------------------------------------------------------------------
	if (cl->energyDep() > 5.) {
	  // poor-man's translation
	  //	  vane_id = cl->vaneId();
	  xl      = cl->cog3Vector().x()+3904.1;
	  yl      = cl->cog3Vector().y();

	  e = new TEllipse(xl,yl,50.*cl->energyDep()/100.);
	  e->SetFillStyle(3001);
	  if( geom->hasElement<mu2e::VaneCalorimeter>() ){
	    color = 2;
	  }else if(geom->hasElement<mu2e::DiskCalorimeter>()){
	    color = module_color[cl->vaneId()];
	  }
	  e->SetFillColor(color);
	  e->SetLineColor(color);
	  
	  list_of_ellipses.Add(e);

	  e->Draw();
	}
      }
    }

    // Loop over all straw hits.
    n_displayed_hits = 0;
    for ( size_t ihit=0; ihit<hits.size(); ++ihit ) {
      
      // Data and MC truth for this hit.
      StrawHit        const&   hit(hits.at(ihit));
      //StrawHitMCTruth const& truth(hitsTruth.at(ihit));
      PtrStepPointMCVector  const& mcptr(hits_mcptr.at(ihit));

      _hEnergyDep->Fill( hit.energyDep()/CLHEP::keV );

      // Skip hits with too little energy deposited in the straw.
      if ( hit.energyDep() < minEnergyDep_ ){
        continue;
      }

      // Get the straw information:
      const Straw&             straw = ttracker->getStraw( hit.strawIndex() );
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();

      bool isFromConversion(false);
      for ( size_t j=0; j<mcptr.size(); ++j ){
        StepPointMC const& step = *mcptr.at(j);
	art::Ptr<SimParticle> const& simptr = step.simParticle();
	SimParticleCollection::key_type trackId(step.trackId());
        SimParticle const& sim  = *simptr;
        if ( sim.fromGenerator() ){
          GenParticle* gen = (GenParticle*) &(*sim.genParticle());
          if ( gen->generatorId() == GenId::conversionGun ){
            isFromConversion = true;
            break;
          }
        }
      }
      _hDeltaT->Fill( hit.dt() );
      _hTime->Fill( hit.time() );

      Hep3Vector pos( tt->positionAtZ( mid.z() ) );
      Hep3Vector mom( tt->momentumAtZ( mid.z() ) );
      TwoLinePCA pca( pos, mom.unit(), mid, w );
      //      double dca = pca.dca();

      Hep3Vector posMid( tmid->positionAtZ( mid.z() ) );
      Hep3Vector momMid( tmid->momentumAtZ( mid.z() ) );
      TwoLinePCA pcaMid( posMid, momMid.unit(), mid, w);

      /*
        double dcaMid = pcaMid.dca();
        cout << "Compare: "
        << mid.z() << " "
        << dca << " "
        << dcaMid << " | "
        << dca-dcaMid << " | "
        << mom.unit().dot(w) << " "
        << momMid.unit().dot(w) << " | "
        << pos << " "
        << posMid <<  " "
        << endl;
      */

      //       if ( isFromConversion ) {
      //         float ntHit[_ntHit->GetNvar()];
      //         _hMissDist->Fill( dca );
      //         ntHit[0] = dca;
      //         ntHit[1] = mid.z();
      //         _ntHit->Fill(ntHit);
      //       }

      // Position along wire, from delta t.
      double v     = trackCal->TimeDiffToDistance( straw.index(), hit.dt() );
      double vnorm = v/straw.getHalfLength();
      double sigv = trackCal->TimeDivisionResolution( straw.index(), vnorm )/2.; // P.Murat

      _hx->Fill(v);
      _hxnorm->Fill(vnorm);

      CLHEP::Hep3Vector x0 = mid + v*w;
      CLHEP::Hep3Vector x1 = x0 + sigv*w;
      CLHEP::Hep3Vector x2 = x0 - sigv*w;

      int color;

      int intime = fabs(hit.time()-event_time) < timeWindow_;

      if ( doDisplay_  ){
        if ( isFromConversion ) {
	  if (intime) color = kRed;
	  else                                   color = kBlue;
	}
	else                                     color = kBlack;
	
	if ((isFromConversion) || (intime)) {
	  line->SetLineColor(color);
	  line->DrawLine( x1.x(), x1.y(), x2.x(), x2.y() );

	  line->DrawLine(x0.x()+5.*w.y(),x0.y()-5*w.x(),x0.x()-5*w.y(),x0.y()+5*w.x());

	  n_displayed_hits++;
	}
      }
    }

    if (doDisplay_) {
      //-----------------------------------------------------------------------------
      // Draw the generated hits.
      //-----------------------------------------------------------------------------
      if (xStep.size() <= 0) {
      }
      else {
	graph = new TGraph( xStep.size(), &xStep[0], &yStep[0]);
	graph->SetMarkerStyle(kOpenTriangleUp);
	graph->SetMarkerSize(0.5);
	graph->Draw("PSAME");
      }

      //-----------------------------------------------------------------------------
      // red marker and arrow: first point on the track.
      //-----------------------------------------------------------------------------
      double xf1 = pos1.x();
      double yf1 = pos1.y();
      TGraph genPoint( 1, &xf1, &yf1 );
      genPoint.SetMarkerColor(kRed);
      genPoint.SetMarkerSize(1.0);
      genPoint.SetMarkerStyle(kFullCircle);
      genPoint.Draw("PSAME");
      
      double arrowLength(200.);
      double xf2 = xf1 + arrowLength*mom1.x()/mom1.perp();
      double yf2 = yf1 + arrowLength*mom1.y()/mom1.perp();
      arrow->SetLineColor(kRed);
      arrow->DrawArrow( xf1, yf1, xf2, yf2, 0.01, ">");
      
      double d0x  = tt->d0x();
      double d0y  = tt->d0y();
      double d0x2 =  tt->d0x() + arrowLength*tt->u0();
      double d0y2 =  tt->d0y() + arrowLength*tt->v0();
      
      arrow->SetLineColor(kBlue);
      arrow->DrawArrow(d0x, d0y, d0x2, d0y2, 0.01, ">");
      //-----------------------------------------------------------------------------
      // blue cross - origin.
      //-----------------------------------------------------------------------------
      double xo = 0.;
      double yo = 0.;
      TGraph genPointo( 1, &xo, &yo );
      genPointo.SetMarkerColor(kBlue);
      genPointo.SetMarkerSize(2.5);
      genPointo.SetMarkerStyle(kPlus);
      genPointo.Draw("PSAME");
      
      canvas_->Modified();
      canvas_->Update();
      
      printf("N(hits) = %5i, N(displayed hits): %5i\n",(int) hits.size(),n_displayed_hits);
      printf("number of clusters      : %5i\n",fNClusters);
      
      char junk(0);
      
      if ( clickToAdvance_ ) {
        cerr << "Double click in the Canvas " << moduleLabel_ << " to continue:" ;
        gPad->WaitPrimitive();
      }
      else{
        cerr << "Enter any character to continue: ";
        cin >> junk;
      }
      cerr << endl;
    }
    //    END_OF_ROUTINE:;
    //-----------------------------------------------------------------------------
    // memory cleanup
    //-----------------------------------------------------------------------------
    list_of_ellipses.Delete();
    
    if (tt) {
      delete tt;
      delete tg;
      delete tmid;
    }

    if (graph) delete graph;
    
  } // end of ::analyze.
    
}

using mu2e::HitDisplay;
DEFINE_ART_MODULE(HitDisplay);
