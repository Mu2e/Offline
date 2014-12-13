///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: MuHitDisplay_module.cc,v 1.6 2014/09/20 17:54:06 murat Exp $
// $Author: murat $
// $Date: 2014/09/20 17:54:06 $
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"
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
#include "KalmanTests/inc/TrkFitDirection.hh"
#include "TrkBase/TrkParticle.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepPtrCollection.hh"

#include "CalPatRec/inc/CalTimePeak.hh"

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
#include "TDatabasePDG.h"

// Other includes
// #include "CLHEP/Units/SystemOfUnits.h"

#include "Stntuple/gui/gui/TStnVisManager.hh"
					// nodes represent objects to be displayed
#include "Stntuple/gui/THeaderVisNode.hh"
#include "Stntuple/gui/TCalVisNode.hh"
#include "Stntuple/gui/TTrkVisNode.hh"
#include "Stntuple/gui/TStrawHitVisNode.hh"
#include "Stntuple/gui/TMcTruthVisNode.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "Stntuple/mod/TAnaDump.hh"
#include "Stntuple/mod/THistModule.hh"

#include "Stntuple/obj/TStnHeaderBlock.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class MuHitDisplay : public THistModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string        moduleLabel_;	             // this module label
    std::string        processName_;
    std::string        generatorModuleLabel_; 
    std::string        fG4ModuleLabel;
    std::string        caloClusterModuleLabel_;
    
    std::string        fCaloClusterAlgorithm;
    std::string        fCaloClusterSeeding;
    std::string        producerName_;
    std::string        fStrawHitMaker;
    std::string        fStrawHitPosMaker;
    std::string        fStrawHitFlagMaker;
    std::string        fTrkPatRecModuleLabel;
    std::string        fCrystalHitMaker;
    std::string        fTrkExtrapol;
    std::string        fTrkCalMatch;
    std::string        fPidModuleLabel;

    TrkFitDirection    fTrkDirection;
    TrkParticle        fParticleHypo;


    int                fGeneratorID;
					       // Name of the tracker StepPoint collection
    std::string        trackerStepPoints_;
    
    // Cuts used inside SimParticleWithHits:
    //  - drop hits with too little energy deposited.
    //  - drop SimParticles with too few hits.


    // name for creating module label
    string             fDirectionAndParticle;

    double             minEnergyDep_;
    double             timeWindow_;
    size_t             minHits_;
					// Options to control the display

    bool               fDisplayBackgroundHits;
    bool               clickToAdvance_;
    bool               fPrintHits;
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
    const mu2e::CaloCrystalHitCollection*       fListOfCrystalHits;//
    const mu2e::CaloClusterCollection*          fListOfClusters;   //
    const mu2e::PtrStepPointMCVectorCollection* hits_mcptr;        //
    const mu2e::StepPointMCCollection*          fSteps;            //
    const mu2e::CalTimePeakCollection*          fCalTimePeakColl;  //

    mu2e::SimParticlesWithHits*                 fSimParticlesWithHits; // 

    int       fNClusters;
					// 4 hypotheses: dem, uep, dmm, ump
    int       fNTracks[4];

    const mu2e::KalRepPtrCollection*  fDem; 

    TMarker*            fMarker;

    TStnVisManager*     fVisManager;
    TStnTrackBlock*     fTrackBlock;
    TStnClusterBlock*   fClusterBlock;
    TStnHeaderBlock*    fHeaderBlock;

    TStnTrackID*        fTrackID;

    double              fTrackerP;   	// CE momentum on entrance to the straw tracker
    double              fTrackerPt;	// CE Pt on entrance to the straw tracker

    const CalTimePeak*  fTimePeak;

  public:
    explicit MuHitDisplay(fhicl::ParameterSet const& pset);
    virtual ~MuHitDisplay();

    int      getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
    void     InitVisManager();

    void     printCaloCluster(const CaloCluster* Cl, const char* Opt) ;
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual bool     filter  (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  MuHitDisplay::MuHitDisplay(fhicl::ParameterSet const& pset):
    THistModule(pset,"MuHitDisplay"),
    moduleLabel_              (pset.get<std::string>("module_label"                )),
    processName_(pset.get<string>("processName","")),
    generatorModuleLabel_     (pset.get<std::string>("generatorModuleLabel"        )),
    fG4ModuleLabel            (pset.get<std::string>("g4ModuleLabel"               )),
    caloClusterModuleLabel_   (pset.get<std::string>("caloClusterModuleLabel")),
    fCaloClusterAlgorithm     (pset.get<std::string>("caloClusterAlgorithm"  )),
    fCaloClusterSeeding       (pset.get<std::string>("caloClusterSeeding"    )),

//     producerName_             ("Algo"+mu2e::TOUpper(fCaloClusterAlgorithm)
// 			       +"SeededBy"+mu2e::TOUpper(fCaloClusterSeeding)),
    producerName_             (""),

    fStrawHitMaker            (pset.get<std::string>("strawHitMakerModuleLabel"    )),
    fStrawHitPosMaker         (pset.get<std::string>("strawHitPosMakerModuleLabel" )),
    fStrawHitFlagMaker        (pset.get<std::string>("strawHitFlagMakerModuleLabel")),

    fTrkPatRecModuleLabel     (pset.get<std::string>("trkPatRecModuleLabel"        )),
    fCrystalHitMaker          (pset.get<std::string>("caloCrystalHitsMaker"  )),
    fTrkExtrapol              (pset.get<std::string> ("trkExtrapol"          )),
    fTrkCalMatch              (pset.get<std::string> ("trkCalMatch"          )),
    fPidModuleLabel           (pset.get<std::string> ("pidModuleLabel"       )),

    fTrkDirection             ( (TrkFitDirection::FitDirection)(pset.get<int>("fitDirection"))),
    fParticleHypo             ( (TrkParticle::type)            (pset.get<int>("fitParticle" ))),          

    fGeneratorID              (pset.get<int>        ("generatorID"           ,GenId::conversionGun)),
    trackerStepPoints_        (pset.get<std::string>("trackerStepPoints"           )),
    minEnergyDep_             (pset.get<double>     ("minEnergyDep"         ,0     )),
    timeWindow_               (pset.get<double>     ("timeWindow"           ,1.e6  )),
    minHits_                  (pset.get<unsigned>   ("minHits"                     )),
    fDisplayBackgroundHits    (pset.get<bool>       ("displayBackgroundHits",false )),
    clickToAdvance_           (pset.get<bool>       ("clickToAdvance"       ,false )),
    fPrintHits                (pset.get<bool>       ("printHits"            ,false )),
    fUseStereoHits            (pset.get<bool>       ("useStereoHits"        ,false )),
    fGoodHitMask              (pset.get<std::vector<std::string> >("goodHitMask"   )),
    fBadHitMask               (pset.get<std::vector<std::string> >("badHitMask"    ))
  {

    fApplication          = 0;
    fCanvas               = 0;
    fCal                  = 0;

    fTrackBlock           = new TStnTrackBlock  ();
    fClusterBlock         = new TStnClusterBlock();
    fHeaderBlock          = new TStnHeaderBlock ();

    fVisManager           = TStnVisManager::Instance();
    fSimParticlesWithHits = NULL;

    fTrackID      = new TStnTrackID();

    fDirectionAndParticle = fTrkDirection.name() + fParticleHypo.name(); 
 }

//-----------------------------------------------------------------------------
  MuHitDisplay::~MuHitDisplay() { 
    if (fApplication) delete fApplication; 
  }


//-----------------------------------------------------------------------------
  void MuHitDisplay::beginJob() {

    //     const char oname[] = "MuHitDisplay::beginJob";
    int    tmp_argc(0);
    char** tmp_argv(0);

    if (!gApplication) {
      fApplication = new TApplication( "MuHitDisplay_module", &tmp_argc, tmp_argv );
    }

    // Create a canvas with a guaranteed unique name; the module label is unique within a job.
    TString name  = "canvas_"     + moduleLabel_;
    TString title = "Canvas for " + moduleLabel_;

    fCanvas = new TCanvas(name,title,800,800);
    fMarker = new TMarker(0,0,20);
    fMarker->SetMarkerSize(0.3);

//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------

    const char* charDirectionAndParticle = fDirectionAndParticle.c_str();


    fClusterBlock->AddCollName("mu2e::CaloClusterCollection",
			       caloClusterModuleLabel_.data(),
			       "");
    fClusterBlock->AddCollName("mu2e::TrackClusterLink",fTrkCalMatch.data(),"");
    
    fTrackBlock->AddCollName("mu2e::CaloClusterCollection"         ,caloClusterModuleLabel_.data(),"");
    fTrackBlock->AddCollName("mu2e::KalRepCollection"              ,fTrkPatRecModuleLabel.data()  , charDirectionAndParticle);
    fTrackBlock->AddCollName("mu2e::TrkToCaloExtrapolCollection"   ,fTrkExtrapol.data()           ,"");
    fTrackBlock->AddCollName("mu2e::TrackClusterLink"              ,fTrkCalMatch.data()           ,"");
    fTrackBlock->AddCollName("mu2e::StrawHitCollection"            ,fStrawHitMaker.data()         ,"");
    fTrackBlock->AddCollName("mu2e::PtrStepPointMCVectorCollection",fStrawHitMaker.data()         ,"StrawHitMCPtr");
    fTrackBlock->AddCollName("mu2e::PIDProductCollection"          ,fPidModuleLabel.data()        ,"");

    TAnaDump::Instance()->AddObject("MuHitDisplay::TrackBlock"  ,fTrackBlock  );
    TAnaDump::Instance()->AddObject("MuHitDisplay::ClusterBlock",fClusterBlock);
  }

//-----------------------------------------------------------------------------
// initialize the visualization manager
//-----------------------------------------------------------------------------
  void MuHitDisplay::InitVisManager() {
    const char oname[] = "MuHitDisplay::InitVisManager";

    fVisManager->SetTitleNode(new THeaderVisNode("HeaderVisNode",fHeaderBlock));

    TCalVisNode      *cal_node[2];
    TVisNode         *trk_node;
    TStrawHitVisNode *strh_node;
    TMcTruthVisNode  *mctr_node;
    const mu2e::DiskCalorimeter* dc;

    art::ServiceHandle<mu2e::GeometryService> geom;

    if (geom->hasElement<mu2e::DiskCalorimeter>()) {
      mu2e::GeomHandle<mu2e::DiskCalorimeter> dc_handle;
      dc = dc_handle.get();

      cal_node[0] = new TCalVisNode("CalVisNode#0",&dc->disk(0),0);
      cal_node[0]->SetListOfClusters(&fListOfClusters);
      cal_node[0]->SetListOfCrystalHits(&fListOfCrystalHits);
      cal_node[0]->SetCalTimePeakColl (&fCalTimePeakColl );
      fVisManager->AddNode(cal_node[0]);
      
      cal_node[1] = new TCalVisNode("CalVisNode#1",&dc->disk(1),1);
      cal_node[1]->SetListOfClusters(&fListOfClusters);
      cal_node[1]->SetListOfCrystalHits(&fListOfCrystalHits);
      cal_node[1]->SetCalTimePeakColl (&fCalTimePeakColl );
      fVisManager->AddNode(cal_node[1]);
    }
    else {
      printf("[%s] >>> ERROR : DiskCalorimeter is not defined\n",oname);
      dc = 0;
    }

    trk_node = new TTrkVisNode("TrkVisNode");
    fVisManager->AddNode(trk_node);

    strh_node = new TStrawHitVisNode("StrawHitVisNode");
    strh_node->SetStrawHitColl    (&fStrawHitColl     );
    strh_node->SetStrawHitPosColl (&fStrawHitPosColl  );
    strh_node->SetStrawHitFlagColl(&fStrawHitFlagColl );
    strh_node->SetCalTimePeakColl (&fCalTimePeakColl  );
    strh_node->SetMcPtrColl(&hits_mcptr);
    fVisManager->AddNode(strh_node);

    mctr_node = new TMcTruthVisNode("McTruthVisNode");
    mctr_node->SetListOfHitsMcPtr(&hits_mcptr);
    mctr_node->SetStepPointMCCollection(&fSteps);
    mctr_node->SetSimParticlesWithHits(&fSimParticlesWithHits);
    mctr_node->SetGenpColl(&fGenpColl);
    fVisManager->AddNode(mctr_node);
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  int MuHitDisplay::getData(const art::Event* Evt) {
    const char* oname = "MuHitDisplay::getData";

//-----------------------------------------------------------------------------
//  MC truth - gen particles
//-----------------------------------------------------------------------------
    art::Handle<GenParticleCollection> gensHandle;
    Evt->getByLabel(generatorModuleLabel_,gensHandle);
    if (gensHandle.isValid()) fGenpColl = gensHandle.product();
    else {
      printf(">>> %s ERROR: GenParticleCollection by %s is missing. BAIL OUT\n",
	     oname,generatorModuleLabel_.data());
      return -1;
    }

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    Evt->getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandle);
    if(mcptrHandle.isValid()) hits_mcptr = mcptrHandle.product();
    else {
      printf(">>> %s ERROR:PtrStepPointMCVectorCollection  by %s is missing. BAIL OUT\n",
	     oname,fStrawHitMaker.data());
      return -1;
    }

    art::Handle<StepPointMCCollection> stepsHandle;
    art::Selector getTrackerSteps(art::ProductInstanceNameSelector(trackerStepPoints_) &&
				  art::ProcessNameSelector(processName_) &&
				  art::ModuleLabelSelector(fG4ModuleLabel)  );
    //Evt->getByLabel(fG4ModuleLabel,trackerStepPoints_,stepsHandle);
    Evt->get(getTrackerSteps, stepsHandle);

    if (stepsHandle.isValid()) fSteps = stepsHandle.product();
    else                       fSteps = NULL;
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shH;
    Evt->getByLabel(fStrawHitMaker,shH);
    if (shH.isValid()) fStrawHitColl = shH.product();
    else {
      printf(">>> %s ERROR:StrawHitCollection by %s is missing. BAIL OUT\n",
	     oname,fStrawHitMaker.data());
      return -1;
    }
    
    art::Handle<mu2e::StrawHitPositionCollection> shpH;
    Evt->getByLabel(fStrawHitPosMaker,shpH);
    if (shpH.isValid()) fStrawHitPosColl = shpH.product();
    else {
      printf(">>> %s ERROR:StrawHitPositionCollection by %s is missing. BAIL OUT\n",
	     oname,fStrawHitPosMaker.data());
      return -1;
    }

    art::Handle<mu2e::StrawHitFlagCollection> shfH;
    Evt->getByLabel(fStrawHitFlagMaker,shfH);
    if (shfH.isValid()) fStrawHitFlagColl = shfH.product();
    else {
      printf(">>> %s ERROR: StrawHitFlagCollection by %s is missing. BAIL OUT\n",
	     oname,fStrawHitFlagMaker.data());
      return -1;
    }
//-----------------------------------------------------------------------------
// calorimeter crystal hit data
//-----------------------------------------------------------------------------
    art::Handle<CaloCrystalHitCollection> ccHandle;
    Evt->getByLabel(fCrystalHitMaker.data(),ccHandle);
    if (ccHandle.isValid()) fListOfCrystalHits = (CaloCrystalHitCollection*) ccHandle.product();
    else {
      printf(">>> %s ERROR: CaloCrystalHitCollection by %s is missing. BAIL OUT\n",
	     oname,fCrystalHitMaker.data());
      return -1;
    }
//-----------------------------------------------------------------------------
// calorimeter cluster data
//-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection> calo_cluster_handle;
    Evt->getByLabel(caloClusterModuleLabel_, producerName_ ,calo_cluster_handle);
    if (calo_cluster_handle.isValid()) fListOfClusters = calo_cluster_handle.product();
    else {
      printf(">>> %s ERROR: CaloClusterCollection by %s is missing. BAIL OUT\n",
	     oname,caloClusterModuleLabel_.data());
      return -1;
    }
//-----------------------------------------------------------------------------
// timepeaks 
//-----------------------------------------------------------------------------
    fCalTimePeakColl = NULL;
    fTimePeak        = NULL;

    art::Handle<CalTimePeakCollection> tpch;
    const char* charDirectionAndParticle = fDirectionAndParticle.c_str();
    //   Evt->getByLabel("CalPatRec","DownstreameMinus",tpch);
    Evt->getByLabel("CalPatRec",charDirectionAndParticle,tpch);
    if (tpch.isValid()) { 
      fCalTimePeakColl = tpch.product();
//-----------------------------------------------------------------------------
// find the right time peak to display - display the first one with the track
//-----------------------------------------------------------------------------
      const CalTimePeak* tp;
      int ipeak = -1;
      if (fCalTimePeakColl != NULL) {
	int ntp = fCalTimePeakColl->size();
	for (int i=0; i<ntp; i++) {
	  tp = &fCalTimePeakColl->at(i);
	  if (tp->CprIndex() >= 0) {
	    fTimePeak = tp;
	    ipeak     = i;
	    break;
	  }
	}
      }
      fVisManager->SetTimePeak(ipeak);
    }
//-----------------------------------------------------------------------------
// tracking data - downstream moving electrons
//-----------------------------------------------------------------------------
    art::Handle<KalRepPtrCollection> demHandle;
    Evt->getByLabel(fTrkPatRecModuleLabel.data(),charDirectionAndParticle, demHandle);

    fNTracks[0] = 0;
    fDem        = NULL;
    if (demHandle.isValid()) { 
      fDem        = demHandle.product();
      fNTracks[0] = fDem->size();
    }
//-----------------------------------------------------------------------------
// three other track collections
//-----------------------------------------------------------------------------
//     art::Handle<KalRepPtrCollection> uepHandle;
//     Evt->getByLabel("TrkPatRec2","UpstreamePlus", uepHandle);
//     const KalRepPtrCollection*  uep = &(*uepHandle);
//     fNTracks[1] = uep->size();

    //     art::Handle<KalRepPtrCollection> dmmHandle;
    //     Evt->getByLabel("trkPatRec3","DownstreammuMinus", dmmHandle);
    //     const KalRepPtrCollection*  dmm = &(*dmmHandle);
    //     fNTracks[2] = dmm->size();

    //     art::Handle<KalRepPtrCollection> umpHandle;
    //     Evt->getByLabel("trkPatRec4","UpstreammuPlus", umpHandle);
    //     const KalRepPtrCollection*  ump = &(*umpHandle);
    //     fNTracks[3] = ump->size();

    return 0;
  }

//-----------------------------------------------------------------------------
  void MuHitDisplay::Init(art::Event* Evt) {
//    TStnCluster*    cluster;
    int             id_word, ntrk;
    TStnTrack*      track;

//-----------------------------------------------------------------------------
// initialize tracks and determine track quality
//-----------------------------------------------------------------------------
    StntupleInitMu2eTrackBlock(fTrackBlock,Evt,0);

    fNTracks[0] = fTrackBlock->NTracks();

    ntrk = fNTracks[0];

    for (int i=0; i<ntrk; i++) {
      track          = fTrackBlock->Track(i);
      id_word        = fTrackID->IDWord(track);
      track->fIDWord = id_word;
    }
//-----------------------------------------------------------------------------
// initialize clusters
//-----------------------------------------------------------------------------
    StntupleInitMu2eClusterBlock(fClusterBlock,Evt,0);
    fNClusters = fClusterBlock->NClusters();
//-----------------------------------------------------------------------------
// find momentum of the G4 particle in front of the tracker
//-----------------------------------------------------------------------------
    art::Handle<mu2e::StepPointMCCollection> handle;
    const mu2e::StepPointMCCollection*       coll(0);

    art::Selector  selector(art::ProductInstanceNameSelector("virtualdetector") &&
			    art::ProcessNameSelector        (""         ) && 
			    art::ModuleLabelSelector        ("g4run"    )    );
    Evt->get(selector, handle);

    if (handle.isValid()) coll = handle.product();
    else {
      printf(">>> ERROR in MuHitDisplay::Init: failed to locate StepPointMCCollection for virtual detectors");
      printf(". BAIL OUT. \n");
      return;
    }
    
  int nsteps = coll->size();

  const mu2e::StepPointMC* step;

  const CLHEP::Hep3Vector* mom;

  double px, py, pt;

  fTrackerP  = -1.;
  fTrackerPt = -1.;

  for (int i=0; i<nsteps; i++) {
    step = &coll->at(i);
    if (step->volumeId() == 13) {
      mom = &step->momentum();
      fDump->printStepPointMC(step, "banner+data");

      px = mom->x();
      py = mom->x();

      pt = sqrt(px*px+py*py);

      fTrackerP  = mom->mag();
      fTrackerPt = pt;

      break;
    }
  }
//-----------------------------------------------------------------------------
// then - clusters, tracks are supposed to be already initialized
//-----------------------------------------------------------------------------
//     const mu2e::CaloCluster       *cl;
//     unsigned int nm(0);

//     if (fTrkCalMap) nm = fTrkCalMap->size();

//     for (int i=0; i<fNClusters; i++) {
//       cluster               = fClusterBlock->NewCluster();

//       for(size_t i=0; i<nm; i++) {
// 	//	KalRepPtr const& trkPtr = fTrkCalMap->at(i).first->trk();
// 	//	const KalRep *  const &trk = *trkPtr;

// 	cl = &(*(fTrkCalMap->at(i).second));

// 	if (cl == cluster->fCaloCluster) { 
// 	  cluster->fClosestTrack = fTrackBlock->Track(i);
// 	  break;
// 	}
//       }
//     }
  }


//-----------------------------------------------------------------------------
  bool MuHitDisplay::filter(art::Event& Evt) {
    const char* oname = "MuHitDisplay::filter";

    static int          firstCall (1);
    TText               t;
    TEllipse*           e;
    TGraph*             graph (0);

    TrackTool           *tt(0), *tg(0), *tmid(0);
    const GenParticle*  gen_signal;
    int                 pdg_id, ndisks, rc;
 
    TObjArray           list_of_ellipses;
    int                 n_displayed_hits, color(1), intime;
    size_t              nmc; 

    const int           module_color[2] = {kRed, kMagenta};

    vector<double>      xStep, yStep;

    const KalRep*       trk;
    const CaloCluster   *cl, *clj;
    //    int             vane_id;
    double              xl, yl, event_time;
    float               q;

    TString             opt; 
    Hep3Vector          pos1, mom1, pos2, mom2;

    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

    // Geometry of tracker envelope.

    mu2e::GeomHandle<mu2e::TTracker> ttHandle;
    const mu2e::TTracker* tracker = ttHandle.get();

    TubsParams envelope(tracker->getInnerTrackerEnvelopeParams());

    // Tracker calibration object.
    mu2e::ConditionsHandle<mu2e::TrackerCalibrations> trackCal("ignored");

//-----------------------------------------------------------------------------
// init VisManager - failed to do it in beginJob - what is the right place for doing it?
//-----------------------------------------------------------------------------
    if (firstCall == 1) {
      firstCall = 0;
      InitVisManager();
    }
    fDump->SetEvent(Evt);
//-----------------------------------------------------------------------------
// get event data and initialize data blocks
//-----------------------------------------------------------------------------
    rc = getData(&Evt);

    if (rc < 0) {
      printf(" %s ERROR: not all data products present, BAIL OUT\n",oname);
      return true;
    }

    Init   (&Evt);

    art::ServiceHandle<mu2e::GeometryService> geom;

    // Construct an object that ties together all of the simulated particle and hit info.

    if (fSimParticlesWithHits) delete fSimParticlesWithHits;

    fSimParticlesWithHits = new SimParticlesWithHits( Evt,
						      fG4ModuleLabel,
						      fStrawHitMaker,
						      trackerStepPoints_,
						      minEnergyDep_,
						      minHits_ );
//-----------------------------------------------------------------------------
// is he finding the first particle?
//-----------------------------------------------------------------------------
    SimParticleCollection::key_type key(1);

    SimParticleInfo const *info(fSimParticlesWithHits->findOrNull(key));
    const StepPointMC     *firstStep; 

    if ( !info ) {
      printf("[%s] *** ERROR: SimParticleInfo missing, SKIPPING EVENT: %10i\n",oname,Evt.event());
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

    //  DISPLAY:;
    
    fCanvas->cd(0);
    fCanvas->Clear();
//-----------------------------------------------------------------------------
// Draw the frame
//-----------------------------------------------------------------------------
    double plotLimits(850.);
    fCanvas->DrawFrame(-plotLimits,-plotLimits,plotLimits,plotLimits);
    
    t.SetText(-800.,900.,Form("[%s] RUN: %10i EVENT: %10i NTRACKS: %4i NCLUSTERS: %4i",
			      oname, Evt.run(),Evt.event(),fNTracks[0],fNClusters));
    t.SetTextSize(0.02);
    t.Draw();

					// Draw the inner and outer arc of the tracker.
    arc->SetLineColor(kBlack);
    arc->DrawArc(0.,0., envelope.outerRadius());
    arc->DrawArc(0.,0., envelope.innerRadius());

    if (fSteps != NULL) {

      SortedStepPoints sortedSteps(key,*fSteps);
    
      const StepPointMC* midStep = & sortedSteps.middleByZ();
    
      if ((midStep == 0) || (firstStep ==0)) {
	printf("[%s] **** ERROR : firstStep = %8p midstep = %8p, BAIL OUT\n",
	       oname,firstStep, midStep);
      }

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
//-----------------------------------------------------------------------------
// The generated signal particle.
//-----------------------------------------------------------------------------
					// this is not necessarily correct (for mixed events) !
      gen_signal = &fGenpColl->at(0);
      pdg_id     = gen_signal->pdgId();
					// ROOT returns charge in units of dbar-quark charge
 
      q    = TDatabasePDG::Instance()->GetParticle(pdg_id)->Charge()/3.;
      tt   = new TrackTool(pdg_id, q,pos1,mom1,1.,Hep3Vector());
      tg   = new TrackTool(pdg_id, q,gen_signal->position(),gen_signal->momentum(),1.,Hep3Vector());
      tmid = new TrackTool(pdg_id, q,pos2  ,mom2 ,1.,Hep3Vector());

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
     
      arc->SetLineColor(kRed);
      arc->DrawArc( tt->xc(), tt->yc(), tt->rho());
      arc->SetLineColor(kMagenta);
      arc->DrawArc( tmid->xc(), tmid->yc(), tmid->rho());
      arc->SetLineColor(kRed);
    }
    double rmin, rmax;
//-----------------------------------------------------------------------------
// Draw disks: the inner and outer calorimeter boundaries
//-----------------------------------------------------------------------------
    mu2e::GeomHandle<mu2e::DiskCalorimeter> dc;
    ndisks = dc->nDisk();

    for (int id=0; id<ndisks; id++) {
      rmin = dc->disk(id).innerRadius();
      rmax = dc->disk(id).outerRadius();
      arccalo->SetLineColor(kGreen+id);
      arccalo->DrawArc(0.,0., rmax);
      arccalo->DrawArc(0.,0., rmin);
    }
//-----------------------------------------------------------------------------
// print reconstructed tracks - all 4 hypotheses for comparison...
//-----------------------------------------------------------------------------
    printf("[%s] NTRACKS = %4i, NCLUSTERS = %4i\n",oname,fNTracks[0],fNClusters);
    printf("\n");
  
    opt = "data";
    if (fPrintHits) opt += "+hits";
    
    if (fNTracks[0] > 0) {
      
      double  d0, om, r, phi0, x0, y0;

      printf(" I  alg ");
      TAnaDump::Instance()->printKalRep(0,"banner");
      for (int i=0; i<fNTracks[0]; i++ ) {
	trk = fDem->at(i).get();
	printf(" %2i dem ",i);
	TAnaDump::Instance()->printKalRep(trk,opt);
//-----------------------------------------------------------------------------
// also display the reconstructed track, use s=0
//-----------------------------------------------------------------------------
//	HelixParams hel = trk->helix(0);
	d0   = trk->helix(0.).d0();
	om   = trk->helix(0.).omega();
	r    = fabs(1./om);
	phi0 = trk->helix(0.).phi0();

	x0   =  -(1/om+d0)*sin(phi0);
	y0   =   (1/om+d0)*cos(phi0);
// 	printf("[MuHitDispla::printHelixParams] d0 = %5.3f r = %5.3f phi0 = %5.3f x0 = %5.3f y0 = %5.3f\n",
// 	       d0, r, phi0, x0, y0);
	
	e = new TEllipse(x0,y0,r);
	e->SetFillStyle(3001);		// make it transparent

	e->SetLineColor(kBlue-7);
	list_of_ellipses.Add(e);
	e->Draw();
      }
    }
//-----------------------------------------------------------------------------
// Calorimeter Clusters: print and prepare to display
//-----------------------------------------------------------------------------
    event_time = -1.;
    if (fNClusters > 0) {
      TAnaDump::Instance()->printCaloCluster(0,"banner");
    
      for (int i=0; i<fNClusters; i++) {
	cl = &fListOfClusters->at(i);
	if (cl->daddy() == -1) {
	  TAnaDump::Instance()->printCaloCluster(cl,opt);
//-----------------------------------------------------------------------------
// event time defined by the first, most energetic cluster
//-----------------------------------------------------------------------------
	  if (i == 0) {
	    event_time = cl->time()+15.;
	  }

	  for (int j=i+1; j<fNClusters; j++) {
	    clj = &fListOfClusters->at(j);
	    if (clj->daddy() == i) {
	      TAnaDump::Instance()->printCaloCluster(clj,opt);
	    }
	  }
	}
      }

      for (int i=0; i<fNClusters; i++) {
	cl = &fListOfClusters->at(i);
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
//-----------------------------------------------------------------------------
// Loop over straw hits. If flagBackgroundHits = true, filter them out
//-----------------------------------------------------------------------------
    int                          ihit, n_straw_hits, display_hit;
    bool                         isFromConversion;
    double                       sigv, vnorm, v, sigr; 
    const StrawHit*              hit; 
    const StrawHitPosition*      hitpos; 
    const StrawHitFlag*          hit_id_word; 
    const CLHEP::Hep3Vector      *mid, *w; 
    const Straw*                 straw; 
    const PtrStepPointMCVector*  mcptr;
    CLHEP::Hep3Vector            vx0, vx1, vx2; 

    n_displayed_hits = 0;

    if (fTimePeak != NULL) n_straw_hits = fTimePeak->NHits();
    else                   n_straw_hits = fStrawHitColl->size();

    for (int ih=0; ih<n_straw_hits; ++ih ) {

      if (fTimePeak != NULL) ihit = fTimePeak->HitIndex(ih);
      else                   ihit = ih;

      hit         = &fStrawHitColl->at(ihit);
      hitpos      = &fStrawHitPosColl->at(ihit);
      hit_id_word = &fStrawHitFlagColl->at(ihit);

      display_hit = 1;

      if (display_hit) {

	mu2e::GenId gen_id(fGeneratorID);
      //StrawHitMCTruth const& truth(hitsTruth.at(ihit));
	mcptr = &hits_mcptr->at(ihit);
	
	// Get the straw information:
	straw = &tracker->getStraw( hit->strawIndex() );
	mid   = &straw->getMidPoint();
	w     = &straw->getDirection();

	isFromConversion = false;

	nmc = mcptr->size();
	for (size_t j=0; j<nmc; ++j ){
	  const art::Ptr<mu2e::StepPointMC>& sptr = mcptr->at(j);
	  const mu2e::StepPointMC* step = sptr.operator ->();

	  //	  SimParticleCollection::key_type trackId(step->trackId());
	  art::Ptr<SimParticle> const& simptr = step->simParticle();
	  const SimParticle* sim  = simptr.operator ->();
	  if (sim == NULL) {
	    printf(">>> ERROR: %s sim == NULL\n",oname);
	  }
	  else {
	    if ( sim->fromGenerator() ){
	      GenParticle* gen = (GenParticle*) &(*sim->genParticle());
	      //	      if ( gen->generatorId() == gen_id /*GenId::conversionGun*/ ){
	      if (gen == gen_signal) {
		isFromConversion = true;
		break;
	      }
	    }
	  }
	}
	
	//	Hep3Vector pos(tt->positionAtZ( mid->z()));
	//	Hep3Vector mom(tt->momentumAtZ( mid->z()));
	//	TwoLinePCA pca(pos, mom.unit(), *mid, *w);
	
	//	Hep3Vector posMid(tmid->positionAtZ( mid->z() ) );
	//	Hep3Vector momMid(tmid->momentumAtZ( mid->z() ) );
	//	TwoLinePCA pcaMid(posMid, momMid.unit(), *mid, *w);
	
	// Position along wire, from delta t.
	
	v     = trackCal->TimeDiffToDistance( straw->index(), hit->dt() );
	vnorm = v/straw->getHalfLength();
	if (fUseStereoHits) {
//-----------------------------------------------------------------------------
// new default, hit position errors come from StrawHitPositionCollection
//-----------------------------------------------------------------------------
	  sigv  = hitpos->posRes(StrawHitPosition::phi); 
	  sigr  = hitpos->posRes(StrawHitPosition::rho); 
	}
	else {
//-----------------------------------------------------------------------------
// old default, draw semi-random errors
//-----------------------------------------------------------------------------
	  sigv  = trackCal->TimeDivisionResolution( straw->index(), vnorm )/2.; // P.Murat
	  sigr  = 5.; // in mm
	}
	
	//	vx0 = (*mid) + v   *(*w);

	vx0 = hitpos->pos();

	vx1 = vx0    + sigv*(*w);
	vx2 = vx0    - sigv*(*w);
	
	intime = fabs(hit->time()-event_time) < timeWindow_;
	
//-----------------------------------------------------------------------------
// choose the hit color
//-----------------------------------------------------------------------------
	color = kBlack;

	if ( isFromConversion ) {
	  if (intime) color = kRed;
	  else        color = kBlue;
	}
	else {
	  //	  if (! hit_id_word->hasAllProperties(fGoodHitMask)) display_hit = 0;
	  if (hit_id_word->hasAnyProperty(fBadHitMask)) {
	    color = kCyan+3;
	  }
	}
	
	if ((isFromConversion) || (intime)) {
// 	  double x,y;
// 	  x =hitpos->pos().x(); 
// 	  y =hitpos->pos().y(); 

// 	  if ((x > 350) && (x < 500) && (y > 0) && (y < 200)) {
// 	    printf(" [MuHitDisplay::filter] hit x,y,z = %10.3f %10.3f %10.3f\n",
// 		   hitpos->pos().x(),hitpos->pos().y(),hitpos->pos().z());
// 	  }

	  line->SetLineColor(color);
	  line->DrawLine(vx1.x(),vx1.y(),vx2.x(),vx2.y() );
	  
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

    if (fSteps != NULL) {
      double d0x  = tt->d0x();
      double d0y  = tt->d0y();
      double d0x2 =  tt->d0x() + arrowLength*tt->u0();
      double d0y2 =  tt->d0y() + arrowLength*tt->v0();
      
      arrow->SetLineColor(kBlue);
      arrow->DrawArrow(d0x, d0y, d0x2, d0y2, 0.01, ">");
    }
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
      
    fVisManager->SetEvent(Evt);
    fVisManager->DisplayEvent();
					// go into interactive mode, till '.q'
    TModule::filter(Evt);
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

    return true;
  } // end of ::filter.
    
}

using mu2e::MuHitDisplay;
DEFINE_ART_MODULE(MuHitDisplay);
