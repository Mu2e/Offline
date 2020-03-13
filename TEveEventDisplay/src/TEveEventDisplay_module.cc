//Author: SMiddleton 
//Date: Jan 2020
//Purpose: Simple version of module to make TEVe based event displays in Offline environment
//This is the first stages of this development 
// ... libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>

// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
#include <TGMsgBox.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
// ... libEG
#include <TParticle.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveStraightLineSet.h>

//#include <sstream>
#include "fstream"

//TEveEventDisplay Headers:
#include  "TEveEventDisplay/src/dict_classes/NavState.h"
#include  "TEveEventDisplay/src/dict_classes/EvtDisplayUtils.h"

// Mu2e Utilities
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

//Mu2e Tracker Geom:
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrkDiag/inc/ComboHitInfo.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

//Collections:
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
// Mu2e diagnostics
using namespace std;
using namespace mu2e;

void setRecursiveColorTransp(TGeoVolume *vol, Int_t color, Int_t transp)
  {
     if(color>=0)vol->SetLineColor(color);
     if(transp>=0)vol->SetTransparency(transp);
     Int_t nd = vol->GetNdaughters();
     for (Int_t i=0; i<nd; i++) {
        setRecursiveColorTransp(vol->GetNode(i)->GetVolume(), color, transp);
     }
  }

namespace mu2e 
{
  class TEveEventDisplay : public art::EDAnalyzer {
    public:

struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<int> mcdiag{Name("mcdiag"), Comment("set on for MC info"),2};
	      fhicl::Atom<art::InputTag> strawDigisTag{Name("strawDigisTag"),Comment("strawDigiTag")};
	      
   	      fhicl::Atom<std::string> g4ModuleLabel{Name("g4ModuleLabel"), Comment("")};
	      fhicl::Atom<double> minEnergyDep{Name("minEnergyDep"), Comment("choose minium energy"), 50};
	       fhicl::Atom<int> minHits{Name("minHits"), Comment(""), 2};
	       fhicl::Atom<bool> doDisplay{Name("doDisplay"), Comment(""), true};
	       fhicl::Atom<bool> clickToAdvance{Name("clickToAdvance"), Comment(""), true}; 
               fhicl::Atom<bool> showEvent{Name("showEvent"), Comment(""),true};     
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;
    explicit TEveEventDisplay(const Parameters& conf);
     virtual ~TEveEventDisplay();
     virtual void beginJob();
     virtual void beginRun(const art::Run& run);
     void analyze(const art::Event& e);
     virtual void endJob();
 private:
     Config _conf;
     int _mcdiag;
     Int_t _evt; 
    
     const StrawDigiCollection* _stcol;
     art::InputTag strawDigisTag_;
     std::string g4ModuleLabel_;
     
     // Cuts used inside SimParticleWithHits:
     //  - drop hits with too little energy deposited.
     //  - drop SimParticles with too few hits.
      double minEnergyDep_;
      size_t minHits_;

      bool doDisplay_;
      bool clickToAdvance_;
      bool showEvent_;
      TApplication* application_;
      TDirectory*   directory_ = nullptr;
      bool            drawGenTracks_;
      bool            drawHits_;
      Double_t        hitMarkerSize_;
      Double_t        trkMaxR_;
      Double_t        trkMaxZ_;
      Double_t        trkMaxStepSize_;
      Double_t        camRotateCenterH_;
      Double_t        camRotateCenterV_;
      Double_t        camDollyDelta_;


      TGTextEntry      *fTeRun,*fTeEvt;
      TGLabel          *fTlRun,*fTlEvt;
   
      EvtDisplayUtils *visutil_ = new EvtDisplayUtils();
      
      bool foundEvent = false;
      void MakeNavPanel();
      bool FindData(const art::Event& event);
};

TEveEventDisplay::TEveEventDisplay(const Parameters& conf) :
	art::EDAnalyzer(conf),
	_mcdiag	(conf().mcdiag()),
        strawDigisTag_(conf().strawDigisTag()),
        g4ModuleLabel_(conf().g4ModuleLabel()),
        minEnergyDep_(conf().minEnergyDep()),
        minHits_(conf().minHits()),
	doDisplay_(conf().doDisplay()),
        clickToAdvance_(conf().clickToAdvance()),
        showEvent_(conf().showEvent()){
		visutil_ = new EvtDisplayUtils();
	}


TEveEventDisplay::~TEveEventDisplay(){}

/*-------Create Control Panel For Event Navigation----""*/
void TEveEventDisplay::MakeNavPanel()
{
  cout<<"Here the is is usually code to make Event and Run Selection"<<endl;
}

void TEveEventDisplay::beginJob(){
  cout<<"Beginning Job ... this has been deleted to help debugging...."<<endl;

}


void TEveEventDisplay::beginRun(const art::Run& run){
  cout<<"importing GDML "<<endl;
  // Import the GDML of entire Mu2e Geometry
  TGeoManager* geom = new TGeoManager("geom","Geom");
  geom->TGeoManager::SetNavigatorsLock(kFALSE);
  geom->TGeoManager::ClearNavigators();
  cout<<"the following causes the break:geom = geom-> TGeoManager::Import(TEveEventDisplay/src/mu2e.gdml...."<<endl;
  geom = geom->TGeoManager::Import("TEveEventDisplay/src/mu2e.gdml");
 
  //Get Top Volume
  TGeoVolume* topvol = geom->GetTopVolume();
  
  //Set Top Volume for gGeoManager:
  gGeoManager->SetTopVolume(topvol);
  gGeoManager->SetTopVisible(kTRUE);
  int nn = gGeoManager->GetNNodes();
  printf("nodes in geom = %d\n",nn);
  //Get Top Node:
  TGeoNode* topnode = gGeoManager->GetTopNode();
  TEveGeoTopNode* etopnode = new TEveGeoTopNode(gGeoManager, topnode);
  etopnode->SetVisLevel(4);
  etopnode->GetNode()->GetVolume()->SetVisibility(kFALSE);
  //Set colours to allow transparency:
  setRecursiveColorTransp(etopnode->GetNode()->GetVolume(), kWhite-10,70);
  gEve->AddGlobalElement(etopnode);
  
}


void TEveEventDisplay::analyze(const art::Event& event){
 
  _evt = event.id().event();
  cout<<" This has been removed for debugging"<<endl;
  if(showEvent_ ){
  	FindData(event);
  }
  
} 


bool TEveEventDisplay::FindData(const art::Event& evt){
	_stcol = 0; 
        cout<<"Test to see if we can retrieve data..."<<endl;
	auto chH = evt.getValidHandle<mu2e::StrawDigiCollection>(strawDigisTag_);
	_stcol = chH.product();
	foundEvent = true;
	return _stcol != 0;
       }

void TEveEventDisplay::endJob(){
	if(!foundEvent){
		char msg[300];
		sprintf(msg, "Reached end of file but #%i has not been found", true);
	        new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), "Event Not Found", msg, kMBIconExclamation,kMBOk);
	}

}  
	
}

using mu2e::TEveEventDisplay;
DEFINE_ART_MODULE(TEveEventDisplay);
