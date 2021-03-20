// ... libCore
#include <TApplication.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>
#include <TString.h>
#include <TRandom.h>
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
#include <TEveBox.h>
// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include <iostream>
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
  class TEveNonGDMLTest : public art::EDAnalyzer {
	public:

      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};   
      };

      typedef art::EDAnalyzer::Table<Config> Parameters;
      explicit TEveNonGDMLTest(const Parameters& conf);
      virtual ~TEveNonGDMLTest();
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
      TEveBox* b = new TEveBox;
      fhicl::ParameterSet _pset;
      void MakeTEveMu2eMainWindow();
  };

  TEveNonGDMLTest::TEveNonGDMLTest(const Parameters& conf) :
  art::EDAnalyzer(conf),
  _diagLevel(conf().diagLevel())
  
	{}

  TEveNonGDMLTest::~TEveNonGDMLTest(){}

  
  void TEveNonGDMLTest::beginJob(){
    std::cout<<"[Starting TEveNonGDMLTest::beginJob()]"<<std::endl;
    directory_ = gDirectory;
    if ( !gApplication ){
      int    tmp_argc(0);
      char** tmp_argv(0);
      application_ = new TApplication( "noapplication", &tmp_argc, tmp_argv );
    }
    // Initialize global Eve application manager (return gEve)
    TEveManager::Create();
    std::cout<<"Point 1"<<std::endl;
    TRandom& r = * gRandom;
    b->SetMainColor(kCyan);
    b->SetMainTransparency(0);
     std::cout<<"Point 2"<<std::endl;
    Float_t x=0;
    Float_t y=0;
    Float_t z=0;
    double xnew = r.Uniform(-5,5);
    double ynew = r.Uniform(-5,5);
    double znew = r.Uniform(-5,5);
     std::cout<<"Point 3"<<xnew<<" "<<ynew<<" "<<znew<<std::endl;
    b->SetVertex(0, x - 10 + xnew, y - 10 + ynew, z - 10 + znew);
    b->SetVertex(1, x - 10 + xnew, y + 10 + ynew, z - 10 + znew);
    b->SetVertex(2, x + 10 + xnew, y + 10 + ynew, z - 10 + znew);
    b->SetVertex(3, x + 10 + xnew, y - 10 + ynew, z - 10 + znew);
    b->SetVertex(4, x - 10 + xnew, y - 10 + ynew, z + 10 + znew);
    b->SetVertex(5, x - 10 + xnew, y + 10 + ynew, z + 10 + znew);
    b->SetVertex(6, x + 10 + xnew, y + 10 + ynew, z + 10 + znew);
    b->SetVertex(7, x + 10 + xnew, y - 10 + ynew, z + 10 + znew);
    std::cout<<"Point 4"<<std::endl;
    gEve->AddElement(b);
    gEve->Redraw3D(kTRUE);
    std::cout<<"Point 5"<<std::endl;
    std::cout<<"[Ending TEveNonGDMLTest::beginJob()]"<<std::endl;
  }

  void TEveNonGDMLTest::beginRun(const art::Run& run){
    std::cout<<"[Starting TEveNonGDMLTest::beginRun()]"<<std::endl;
    std::cout<<"[Ending TEveNonGDMLTest::beginRun()]"<<std::endl;
  }

  void TEveNonGDMLTest::analyze(const art::Event& event){
    std::cout<<"[In TEveNonGDMLTest::analyze()]"<<std::endl;

    int enter;
    std::cout<<"Press 1 to end test: "<<std::endl;
    std::cin>>enter;
    std::cout<<"[Ending TEveNonGDMLTest::analyze()]"<<std::endl;
  } 


  void TEveNonGDMLTest::endJob(){}  

}
using mu2e::TEveNonGDMLTest;
DEFINE_ART_MODULE(TEveNonGDMLTest);
