// ... libCore
#include <TApplication.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>
#include <TString.h>
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
#include <TGLViewer.h>
#include <TGIcon.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TGeoNode.h>
#include <TGeoPhysicalNode.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

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
  class TEveGDMLTest : public art::EDAnalyzer {
	public:

      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};   
      };

      typedef art::EDAnalyzer::Table<Config> Parameters;
      explicit TEveGDMLTest(const Parameters& conf);
      virtual ~TEveGDMLTest();
      virtual void beginJob() override;
      virtual void beginRun(const art::Run& run) override;
      virtual void analyze(const art::Event& e);
      virtual void endJob() override;
    private:
      Config _conf;
      int _diagLevel;     
      bool isFirstEvent = true;
      TApplication* application_;
      TDirectory*   directory_ = nullptr;   
      TGeoManager *geom;
      fhicl::ParameterSet _pset;
      void MakeTEveMu2eMainWindow();
  };

  TEveGDMLTest::TEveGDMLTest(const Parameters& conf) :
  art::EDAnalyzer(conf),
  _diagLevel(conf().diagLevel())
  
	{}

  TEveGDMLTest::~TEveGDMLTest(){}
  
  void TEveGDMLTest::MakeTEveMu2eMainWindow()
  {
    TEveBrowser* browser = gEve->GetBrowser();
    browser->StartEmbedding(TRootBrowser::kLeft); 
    TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
    frmMain->SetWindowName("EVT NAV");
    frmMain->SetCleanup(kDeepCleanup);
    TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
    {
      TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );
      TGPictureButton* b = 0;
      b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoBack.gif"));
      navFrame->AddFrame(b);

      TGPictureButton* f = 0;
      f = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoForward.gif"));
      navFrame->AddFrame(f);

      std::string logoFile = "TEveEventDisplay/src/Icons/mu2e_logo_oval.png";
      const TGPicture *logo = gClient->GetPicture(logoFile.c_str());
      TGIcon *icon = new TGIcon(navFrame,logo,50,50);
      navFrame->AddFrame(icon,new TGLayoutHints(kLHintsLeft,20,0,0,0));

      frmMain->AddFrame(navFrame);
      TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
      frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
      frmMain->MapSubwindows();
      frmMain->Resize();
      frmMain->MapWindow();
      browser->StopEmbedding();
      browser->SetTabTitle("Event Nav", 0);
    }
  }
  
  void TEveGDMLTest::beginJob(){
    std::cout<<"[Starting TEveGDMLTest::beginJob()]"<<std::endl;
    directory_ = gDirectory;
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
    }
    // Initialize global Eve application manager (return gEve)
    TEveManager::Create();
    MakeTEveMu2eMainWindow();
    gEve->AddEvent(new TEveEventManager("Event", "Empty Event"));
    TGLViewer *glv = gEve->GetDefaultGLViewer();
    glv->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
    std::cout<<"[Ending TEveGDMLTest::beginJob()]"<<std::endl;
  }

  void TEveGDMLTest::beginRun(const art::Run& run){
    std::cout<<"[Starting TEveGDMLTest::beginRun()]"<<std::endl;
    
    if(gGeoManager){
      gGeoManager->GetListOfNodes()->Delete();
      gGeoManager->GetListOfVolumes()->Delete();
      gGeoManager->GetListOfShapes()->Delete();
    }
    gEve->GetGlobalScene()->DestroyElements();
    
    // Import the GDML of entire Mu2e Geometry
    geom = geom->TGeoManager::Import("TEveEventDisplay/src/fix.gdml");

    //Get Top Volume
    TGeoVolume* topvol = geom->GetTopVolume();

    //Set Top Volume for gGeoManager:
    gGeoManager->SetTopVolume(topvol);
    //Set visibility:
    gGeoManager->SetTopVisible(kTRUE);//HERE--> the issue!!

    //Get Top Node:
    TGeoNode* topnode = gGeoManager->GetTopNode();
    TEveGeoTopNode* etopnode = new TEveGeoTopNode(gGeoManager, topnode);
    etopnode->SetVisLevel(4);
    etopnode->GetNode()->GetVolume()->SetVisibility(kFALSE);
    setRecursiveColorTransp(etopnode->GetNode()->GetVolume(), kWhite-10,70);
    std::cout<<"Point 5"<<std::endl;
    //Add static detector geometry to global scene
    gEve->AddGlobalElement(etopnode);
    geom->Draw("ogl");
    std::cout<<"[Ending TEveGDMLTest::beginRun()]"<<std::endl;
  }

  void TEveGDMLTest::analyze(const art::Event& event){
    std::cout<<"[In TEveGDMLTest::analyze()]"<<std::endl;
    int eventid = event.id().event();
    int runid = event.run();
    int subrunid = event.subRun();
    std::cout<<"Drawing Run : "<<runid<<" Sub-Run "<<subrunid<<" Event : "<<eventid<<std::endl;
    
    if(!isFirstEvent){
      gEve->GetViewers()->DeleteAnnotations();
      gEve->GetCurrentEvent()->DestroyElements();
    }
    // Import event into ortho views and apply projections
    //TEveElement* currevt = gEve->GetCurrentEvent();

    geom->Draw("ogl");
    gPad->WaitPrimitive();
    isFirstEvent = false;
    std::cout<<"[Ending TEveGDMLTest::analyze()]"<<std::endl;
  } 


  void TEveGDMLTest::endJob(){}  

}
using mu2e::TEveGDMLTest;
DEFINE_ART_MODULE(TEveGDMLTest);
