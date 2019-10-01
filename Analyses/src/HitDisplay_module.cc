///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display.
//
// $Id: HitDisplay_module.cc,v 1.24 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Contact person:  Pavel Murat, Gianantonio Pezzulo
//
// What this event display shows: 2D view of the straw hits, tracks, and calorimeter clusters
//
// Straw hit display mode:
// -----------------------
// displayBackgroundHits : false - respect hit flags set by FlagStrawHits_module (default)
//                         the timing hit flag check is commented out - I didn't have a MC file
//                         in hands to check
//                       : true  - display all hits
// useStereoHits         : true : displayed hit position errors are defined by the StrawHitPositionCollection
//                       : false: semi-random assignment (sigr=5mm, sigp = strawTimeDivisionErr/2.)
//
// small black triangles: MC truth on the trajectory
//
// red   hits: hits on a track reconstructed as a downstream-moving conversion electron within the time window
//             size of the timing window can be redefined via  talk-to (input .fcl file)
// blue  hits: hits of the conversion electron track, which fall outside the time window
// black hits: hits produced by anything, but the conversion electron
//
// a hit is displayed by a cross with the radial error bar of 5mm and the resolution
// along the straw of sigma(time division)/2, to better guide the eye
//
// green circles - contour of the disk-based calorimeter
// clusters on the 1st disk are shown in red, on the second disk - in pink
//
// there are few other features which we need to document
//
// .fcl file to use: Analyses/test/hitDisplay.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"

#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"

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
//-----------------------------------------------------------------------------
// Module labels
//-----------------------------------------------------------------------------
    std::string        moduleLabel_;	             // this module label
    std::string        generatorModuleLabel_;
    std::string        g4ModuleLabel_;
    std::string        fStrawHitMaker;
    std::string        fStrawHitPosMaker;
    std::string        fStrawHitFlagMaker;

    // Name of the tracker StepPoint collection
    std::string        trackerStepPoints_;

    // Cuts used inside SimParticleWithHits:
    //  - drop hits with too little energy deposited.
    //  - drop SimParticles with too few hits.

    double             minEnergyDep_;
    double             timeWindow_;
    size_t             minHits_;
					// Options to control the display

    bool               fDisplayBackgroundHits;
    bool               clickToAdvance_;
    bool               printHits_;
    int                fUseStereoHits;
					// hit flag bits which should be ON and OFF
    mu2e::StrawHitFlag         fGoodHitMask;
    mu2e::StrawHitFlag         fBadHitMask;

    TApplication*        fApplication;
    TCanvas*             fCanvas;

    const mu2e::Calorimeter*                    fCal;              //

    const mu2e::StrawHitCollection*             fStrawHitColl;     //
    const mu2e::GenParticleCollection*          fGenpColl;         //
    const mu2e::StrawHitPositionCollection*     fStrawHitPosColl;  //
    const mu2e::StrawHitFlagCollection*         fStrawHitFlagColl; //
    const mu2e::CaloClusterCollection*          fListOfClusters;   //
    const mu2e::PtrStepPointMCVectorCollection* hits_mcptr;
    const mu2e::StepPointMCCollection*          fSteps;

    int       fNClusters;
					// 4 hypotheses: dem, uep, dmm, ump
    int       fNTracks[4];

    const mu2e::KalRepCollection*  fDem;

    TMarker*  fMarker;

    // Tracker conditions object.
    ProditionsHandle<StrawResponse> strawResponse_h;

  public:
    explicit HitDisplay(fhicl::ParameterSet const& pset);
    virtual ~HitDisplay();

    void     getData(const art::Event* Evt);

    void     printCaloCluster(const CaloCluster* Cl, const char* Opt) ;
    void     printKalRep(const KalRep* Trk, const char* Opt);
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    void     beginJob( );
    void     analyze(art::Event const& e );
  };


  HitDisplay::HitDisplay(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    moduleLabel_              (pset.get<std::string>("module_label"                )),
    generatorModuleLabel_     (pset.get<std::string>("generatorModuleLabel"        )),
    g4ModuleLabel_            (pset.get<std::string>("g4ModuleLabel"               )),
    fStrawHitMaker            (pset.get<std::string>("comboHitMakerModuleLabel"    )),
    fStrawHitPosMaker         (pset.get<std::string>("comboHitPosMakerModuleLabel" )),
    fStrawHitFlagMaker        (pset.get<std::string>("comboHitFlagMakerModuleLabel")),
    trackerStepPoints_        (pset.get<std::string>("trackerStepPoints"           )),
    minEnergyDep_             (pset.get<double>     ("minEnergyDep",0              )),
    timeWindow_               (pset.get<double>     ("timeWindow"  ,1.e6           )),
    minHits_                  (pset.get<unsigned>   ("minHits"                     )),
    fDisplayBackgroundHits    (pset.get<bool>       ("displayBackgroundHits",false )),
    clickToAdvance_           (pset.get<bool>       ("clickToAdvance"       ,false )),
    printHits_                (pset.get<bool>       ("printHits"            ,false )),
    fUseStereoHits            (pset.get<bool>       ("useStereoHits"        ,false )),
    fGoodHitMask              (pset.get<std::vector<std::string> >("goodHitMask"   )),
    fBadHitMask               (pset.get<std::vector<std::string> >("badHitMask"    ))
  {

    fApplication = 0;
    fCanvas      = 0;
    fCal         = 0;
  }

//-----------------------------------------------------------------------------
  HitDisplay::~HitDisplay() {
    if (fApplication) delete fApplication;
  }


  //-----------------------------------------------------------------------------
  void HitDisplay::beginJob(){
    int    tmp_argc(0);
    char** tmp_argv(0);

    if (!gApplication) {
      fApplication = new TApplication( "HitDisplay_module", &tmp_argc, tmp_argv );
    }

    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + moduleLabel_;
    TString title = "Canvas for " + moduleLabel_;

    fCanvas = new TCanvas(name,title,800,800);
    fMarker = new TMarker(0,0,20);
    fMarker->SetMarkerSize(0.3);
  }


//-----------------------------------------------------------------------------
// 'banner' : print banner
// 'data'   : print track data
// 'hits'   : print hits
//-----------------------------------------------------------------------------
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
	TrkStrawHitVector tshv;
	convert(Trk->hitVector(),tshv);

	printf("---------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");
	printf(" ih U A     len      rms       x          y          z       HitT     HitDt");
	printf("  SInd  Dev Sec Lay  N  Iamb     T0    Rdrift     Xs         Ys          Zs        resid\n");
	printf("---------------------------------------------------------------------------");
	printf("-------------------------------------------------------------------------------\n");

	int i(0);
        for(auto ihit=tshv.begin(); ihit!= tshv.end(); ++ihit){
	  //	  const KalHit* kh = (*it).kalHit();

	  // TrkStrawHit inherits from TrkHitOnTrk

  	  const TrkStrawHit* hit = (*ihit);

	  const ComboHit* sh = &hit->comboHit();
	  Straw*   straw = (Straw*) &hit->straw();

	  double len = hit->fltLen();

	  HepPoint  plen = Trk->position(len);

	  printf("%3i %1i %1i %10.3f %6.3f %10.3f %10.3f %10.3f %8.3f %7.3f",
		 ++i,
		 hit->hitFlag(),
		 hit->isActive(),
		 len,
		 hit->hitRms(),
		 plen.x(),plen.y(),plen.z(),
		 sh->time(), sh->wireDist()
		 );

	  Hep3Vector pos;
	  hit->hitPosition(pos);
	  printf("%6i %3i %3i %3i %3i %3i %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f\n",
		 straw->id().asUint16(),
		 straw->id().getPlane(),
		 straw->id().getPanel(),
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

    typedef art::Ptr< CaloCrystalHit>                CaloCrystalHitPtr;
    typedef std::vector<CaloCrystalHitPtr>     CaloCrystalHitPtrVector;

    const char oname[] = "HitDisplay::printCaloCluster";

    int     ir, iz;
    TString opt = Opt;

    art::ServiceHandle<GeometryService> geom;
    GeomHandle<Calorimeter> cg;

    printf("[%s] >>> ERROR : yet to learn how to print the crystal indices in case of disks\n", oname);

    if ((opt == "") || (opt.Index("banner") >= 0)) {
      printf("-----------------------------------------------------------------------\n");
      printf("DISKID       Time   Row   Col   Energy       X          Y        Z     \n");
      printf("-----------------------------------------------------------------------\n");
    }

    if ((opt == "") || (opt.Index("data") >= 0)) {

      const CaloCrystalHitPtrVector caloClusterHits = Cl->caloCrystalHitsPtrVector();

      int nh = caloClusterHits.size();
      Hep3Vector pos;
      //-----------------------------------------------------------------------------
      // print individual crystals in local disk coordinate system
      //-----------------------------------------------------------------------------
      for (int i=0; i<nh; i++) {
	const CaloCrystalHit* hit = &(*caloClusterHits.at(i));
	int id = hit->id();

	pos = cg->crystal(id).localPosition();


 	iz = -1;
	ir = -1;




	printf("%6i  %10.3f %5i %5i %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	       id,
	       hit->time(),
	       iz,ir,
	       hit->energyDep(),
	       pos.x(),
	       pos.y(),
	       pos.z(),
	       hit->energyDepTot()
	       );
      }
    }
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void HitDisplay::getData(const art::Event* Evt) {
    //    const char* oname = "HitDisplay::getData";

//-----------------------------------------------------------------------------
//  MC truth - gen particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> gensHandle;
    Evt->getByLabel(generatorModuleLabel_,gensHandle);
    fGenpColl = gensHandle.product();

    // art::Handle<SimParticleCollection> simsHandle;
    // Evt->getByLabel(g4ModuleLabel_,simsHandle);
    // SimParticleCollection const& sims = *simsHandle;

    // art::Handle<StrawHitMCTruthCollection> hitsTruthHandle;
    // Evt->getByLabel(hitMakerModuleLabel_,hitsTruthHandle);
    // StrawHitMCTruthCollection const& hitsTruth = *hitsTruthHandle;

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    Evt->getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandle);
    hits_mcptr = mcptrHandle.product();

    art::Handle<StepPointMCCollection> stepsHandle;
    Evt->getByLabel(g4ModuleLabel_,trackerStepPoints_,stepsHandle);
    fSteps = stepsHandle.product();
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shH;
    Evt->getByLabel(fStrawHitMaker,shH);
    fStrawHitColl = shH.product();

    art::Handle<mu2e::StrawHitPositionCollection> shpH;
    Evt->getByLabel(fStrawHitPosMaker,shpH);
    fStrawHitPosColl = shpH.product();

    art::Handle<mu2e::StrawHitFlagCollection> shfH;
    Evt->getByLabel(fStrawHitFlagMaker,shfH);
    fStrawHitFlagColl = shfH.product();
//-----------------------------------------------------------------------------
// calorimeter cluster data
//-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection> calo_cluster_handle;
    Evt->getByLabel("makeCaloCluster","AlgoCLOSESTSeededByENERGY",calo_cluster_handle);
    fListOfClusters = (CaloClusterCollection*) &(*calo_cluster_handle);
    fNClusters      = fListOfClusters->size();
//-----------------------------------------------------------------------------
// tracking data - upstream moving electrons
//-----------------------------------------------------------------------------
    art::Handle<KalRepCollection> demHandle;
    Evt->getByLabel("trkPatRec1","DownstreameMinus", demHandle);
    fDem = demHandle.product();
    fNTracks[0] = fDem->size();
//-----------------------------------------------------------------------------
// three other track collections
//-----------------------------------------------------------------------------
    art::Handle<KalRepCollection> uepHandle;
    Evt->getByLabel("trkPatRec2","UpstreamePlus", uepHandle);
    const KalRepCollection*  uep = &(*uepHandle);
    fNTracks[1] = uep->size();

    //     art::Handle<KalRepCollection> dmmHandle;
    //     Evt->getByLabel("trkPatRec3","DownstreammuMinus", dmmHandle);
    //     const KalRepCollection*  dmm = &(*dmmHandle);
    //     fNTracks[2] = dmm->size();

    //     art::Handle<KalRepCollection> umpHandle;
    //     Evt->getByLabel("trkPatRec4","UpstreammuPlus", umpHandle);
    //     const KalRepCollection*  ump = &(*umpHandle);
    //     fNTracks[3] = ump->size();

  }

  //-----------------------------------------------------------------------------
  void HitDisplay::analyze(art::Event const& Evt) {
    const char* name = "HitDisplay::analyze";

    TText          t;
    TEllipse*      e;
    TGraph*        graph (0);

    TrackTool      *tt(0), *tg(0), *tmid(0);
    const GenParticle* gen;

    TObjArray      list_of_ellipses;
    int            n_displayed_hits, color(-1), intime;
    size_t         nmc;

    const int module_color[2] = {kRed, kMagenta};

    vector<double> xStep,yStep;

    const KalRep*       trk;
    const CaloCluster*  cl;
    double xl, yl,      event_time(-9999.);
    TString             opt;

    printf("[%s] RUN: %10i EVENT: %10i\n",name,Evt.run(),Evt.event());

    // Geometry of tracker envelope.
    GeomHandle<Tracker> tHandle;
    const Tracker* tracker = tHandle.get();

    TubsParams envelope(tracker->getInnerTrackerEnvelopeParams());

    auto const& strawResponse = strawResponse_h.get(Evt.id());

//-----------------------------------------------------------------------------
// get event data
//-----------------------------------------------------------------------------
    getData(&Evt);

    art::ServiceHandle<mu2e::GeometryService> geom;

    // Construct an object that ties together all of the simulated particle and hit info.
    SimParticlesWithHits simsInfo( Evt,
                                   g4ModuleLabel_,
                                   fStrawHitMaker,
                                   trackerStepPoints_,
                                   minEnergyDep_,
                                   minHits_ );

    SimParticleCollection::key_type key(1);

    SimParticleInfo const* info(simsInfo.findOrNull(key));
    const StepPointMC* firstStep;

    if ( !info ) {
      printf("[%s] *** ERROR: SimParticleInfo missing, SKIPPING EVENT: %10i\n",name,Evt.event());
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

    SortedStepPoints sortedSteps(key,*fSteps);

    const StepPointMC* midStep = & sortedSteps.middleByZ();

    if ((midStep == 0) || (firstStep ==0)) {
      printf("[%s] **** ERROR : firstStep = %8p midstep = %8p, BAIL OUT\n",
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
    gen = &fGenpColl->at(0);

    tt   = new TrackTool(gen->pdgId(), -1.,pos1,mom1,1.,Hep3Vector());
    tg   = new TrackTool(gen->pdgId(), -1.,gen->position()      ,gen->momentum()      ,1.,Hep3Vector());
    tmid = new TrackTool(gen->pdgId(), -1.,pos2  ,mom2 ,1.,Hep3Vector());

    int npt = fSteps->size();
    for (int ipt=0; ipt<npt; ++ipt){
      StepPointMC const& step =  fSteps->at(ipt);
      if ( step.totalEDep() > minEnergyDep_ ) {
        xStep.push_back( step.position().x() );
        yStep.push_back( step.position().y() );
      }
    }

    arc->SetFillStyle(0);
    arccalo->SetFillStyle(0);

    //  DISPLAY:;

    fCanvas->cd(0);
    fCanvas->Clear();
//-----------------------------------------------------------------------------
// Draw the frame
//-----------------------------------------------------------------------------
    double plotLimits(850.);
    fCanvas->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);

    t.SetText(-800.,900.,Form("[%s] RUN: %10i EVENT: %10i NTRACKS: %4i NCLUSTERS: %4i",
			      name, Evt.run(),Evt.event(),fNTracks[0],fNClusters));
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
    // draw disks
    //-----------------------------------------------------------------------------
    int nv(0);
    if( geom->hasElement<mu2e::DiskCalorimeter>() ){
      mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
      nv = dc->nDisk();
    }

    double rmin, rmax; //x1(-1.),y1(-1.),x2(-1.),y2(-1.);
    if(geom->hasElement<mu2e::DiskCalorimeter>()) {

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
//-----------------------------------------------------------------------------
// print reconstructed tracks - all 4 hypotheses for comparison...
//-----------------------------------------------------------------------------
    printf("[%s] NTRACKS = %4i, NCLUSTERS = %4i\n",name,fNTracks[0],fNClusters);
    printf("\n");

    opt = "data";
    if (printHits_) opt += "+hits";

    if (fNTracks[0] > 0) {
      printKalRep(0,"banner");
      for (int i=0; i<fNTracks[0]; i++ ) {
	trk = fDem->get(i);
	printf(" %2i dem ",i);
	printKalRep(trk,opt);
      }
    }


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
	  xl      = cl->cog3Vector().x()+3904.1;
	  yl      = cl->cog3Vector().y();

	  e = new TEllipse(xl,yl,50.*cl->energyDep()/100.);
	  e->SetFillStyle(3001);
	  if(geom->hasElement<mu2e::DiskCalorimeter>()){
	    color = module_color[cl->diskId()];
	  }
	  e->SetFillColor(color);
	  e->SetLineColor(color);

	  list_of_ellipses.Add(e);

	  e->Draw();
	}
      }
    }
//-----------------------------------------------------------------------------
// Loop over straw hits. If flagBackgroundHits = true, filter them out
//-----------------------------------------------------------------------------
    int                          n_straw_hits, display_hit;
    bool                         isFromConversion;
    double                       sigv, v, sigr;
    const StrawHit*              hit;
    const StrawHitPosition*      hitpos;
    const StrawHitFlag*          hit_id_word;
    const CLHEP::Hep3Vector      *mid, *w;
    const Straw*                 straw;
    const PtrStepPointMCVector*  mcptr;
    CLHEP::Hep3Vector            vx0, vx1, vx2;

    n_displayed_hits = 0;
    n_straw_hits     = fStrawHitColl->size();

    for (int ihit=0; ihit<n_straw_hits; ++ihit ) {
      hit         = &fStrawHitColl->at(ihit);
      hitpos      = &fStrawHitPosColl->at(ihit);
      hit_id_word = &fStrawHitFlagColl->at(ihit);

      display_hit = 1;
      if (fDisplayBackgroundHits == false) {
	if (! hit_id_word->hasAllProperties(fGoodHitMask)) display_hit = 0;
	if (hit_id_word->hasAnyProperty(fBadHitMask)     ) display_hit = 0;
      }

      if (display_hit) {

      //StrawHitMCTruth const& truth(hitsTruth.at(ihit));
	mcptr = &hits_mcptr->at(ihit);

	// Get the straw information:
	straw = &tracker->getStraw( hit->strawId() );
	mid   = &straw->getMidPoint();
	w     = &straw->getDirection();

	isFromConversion = false;

	nmc = mcptr->size();
	for (size_t j=0; j<nmc; ++j ){
	  const StepPointMC& step = *mcptr->at(j);
	  art::Ptr<SimParticle> const& simptr = step.simParticle();
	  SimParticleCollection::key_type trackId(step.trackId());
	  SimParticle const& sim  = *simptr;
	  if ( sim.fromGenerator() ){
	    GenParticle* gen = (GenParticle*) &(*sim.genParticle());
	    if ( gen->generatorId().isConversion()){
	      isFromConversion = true;
	      break;
	    }
	  }
	}

	Hep3Vector pos(tt->positionAtZ( mid->z()));
	Hep3Vector mom(tt->momentumAtZ( mid->z()));
	TwoLinePCA pca(pos, mom.unit(), *mid, *w);

	Hep3Vector posMid(tmid->positionAtZ( mid->z() ) );
	Hep3Vector momMid(tmid->momentumAtZ( mid->z() ) );
	TwoLinePCA pcaMid(posMid, momMid.unit(), *mid, *w);

	// Position along wire, from delta t.
        double wdist, wderr,halfpv;
        bool vStatus = strawResponse.wireDistance( *straw, hit->energyDep(), hit->dt(), wdist, wderr,halfpv);
        v=wdist;
        sigv=wderr;
        cout << "return status: " << vStatus << endl;

	if (fUseStereoHits) {
//-----------------------------------------------------------------------------
// new default, hit position errors come from StrawHitPositionCollection
//-----------------------------------------------------------------------------
	  sigv  = hitpos->posRes(StrawHitPosition::wire);
	  sigr  = hitpos->posRes(StrawHitPosition::trans);
	}
	else {
//-----------------------------------------------------------------------------
// old default, draw semi-random errors
//-----------------------------------------------------------------------------
	  //sigv is now set in the call to wireDistance
	  sigr  = 5.; // in mm
	}

	vx0 = (*mid) + v   *(*w);
	vx1 = vx0    + sigv*(*w);
	vx2 = vx0    - sigv*(*w);

	intime = fabs(hit->time()-event_time) < timeWindow_;

	if ( isFromConversion ) {
	  if (intime) color = kRed;
	  else        color = kBlue;
	}
	else          color = kBlack;

	if ((isFromConversion) || (intime)) {
	  line->SetLineColor(color);
	  line->DrawLine( vx1.x(), vx1.y(), vx2.x(), vx2.y() );

	  line->DrawLine(vx0.x()+sigr*w->y(),vx0.y()-sigr*w->x(),vx0.x()-sigr*w->y(),vx0.y()+sigr*w->x());
	  n_displayed_hits++;
	}
      }
    }

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

    fCanvas->Modified();
    fCanvas->Update();

    printf("N(hits) = %5i, N(displayed hits): %5i\n",n_straw_hits,n_displayed_hits);
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
