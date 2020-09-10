#include<TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGIcon.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGString.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include<TGLViewer.h>
#include <TGMsgBox.h>
#include<TPolyLine3D.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveParamList.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveStraightLineSet.h>
//TEveMu2e:
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMainWindow.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eHit.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCluster.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCustomHelix.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eCRVEvent.h"
#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eBField.h"

namespace fhicl
{
  class ParameterSet;
}

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

namespace mu2e{

  TEveMu2eMainWindow::TEveMu2eMainWindow() : TGMainFrame(gClient->GetRoot(), 320, 320){}

  TEveMu2eMainWindow::TEveMu2eMainWindow(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet _pset): 
    TGMainFrame(p, w, h),
    fTeRun(0),
    fTeEvt(0),
    fTTEvt(0),
    fTHSlid(0),
    fTlRun(0),
    fTlEvt(0),
    fTlTEvt(0),
    fTlHSlid(0),
    br(0),
    clusterscheck(0),
    hitscheck(0),
    trackscheck(0),
    cosmicscheck(0),
    cosmictrkscheck(0),
    mctrajcheck(0)	
    {
      TEveManager::Create();
      gEve->GetBrowser()->GetTabRight()->SetTab(0);
      gClient->GetRoot();
      TEveBrowser* browser = gEve->GetBrowser();
      FontStruct_t buttonfont = gClient->GetFontByName("-*-helvetica-medium-r-*-*-8-*-*-*-*-*-iso8859-1");
      GCValues_t gval;
      gval.fMask = kGCForeground | kGCFont;
      gval.fFont = gVirtualX->GetFontHandle(buttonfont);
      gClient->GetColorByName("black", gval.fForeground);

      browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

      TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
      frmMain->SetWindowName("EVT NAV");
      frmMain->SetCleanup(kDeepCleanup);

      TGHorizontalFrame* navFrame = new TGHorizontalFrame(frmMain);
      TGVerticalFrame* evtidFrame = new TGVerticalFrame(frmMain);
      {
        TString icondir(TString::Format("%s/icons/", gSystem->Getenv("ROOTSYS")) );

        TGPictureButton *b = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoBack.gif"),1100);
        navFrame->AddFrame(b);
        b->Associate(this);

        TGPictureButton *f = new TGPictureButton(navFrame, gClient->GetPicture(icondir + "GoForward.gif"),1001);
        navFrame->AddFrame(f);
        f->Associate(this);

        TEveParamList *edep = new TEveParamList("CaloEnergySelection");
        gEve->AddToListTree(edep,0);
        edep->AddParameter(TEveParamList::FloatConfig_t("Min Energy Depositied",20,0,100));

        // ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
        fTlRun = new TGLabel(runoFrame,"Run Number");
        fTlRun->SetTextJustify(kTextLeft);
        fTlRun->SetMargins(5,5,5,0);
        runoFrame->AddFrame(fTlRun);

        fTeRun = new TGTextEntry(runoFrame, _runNumber = new TGTextBuffer(5), 1);
        _runNumber->AddText(0, "1");

        runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));
        fTeRun->Associate(this);

        // ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
        fTlEvt = new TGLabel(evnoFrame,"Evt Number");
        fTlEvt->SetTextJustify(kTextLeft);
        fTlEvt->SetMargins(5,5,5,0);
        evnoFrame->AddFrame(fTlEvt);

        fTeEvt = new TGTextEntry(evnoFrame, _eventNumber = new TGTextBuffer(5), 1);
        _eventNumber->AddText(0, "1");

        evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));
        fTeEvt->Associate(this);

        //Create a Time Slider
        TGHorizontalFrame* timeFrame = new TGHorizontalFrame(evtidFrame);
        fTlHSlid = new TGLabel(timeFrame, "Time (ns)");
        fTlHSlid->SetTextJustify(kTextLeft);
        fTlHSlid->SetMargins(5,5,5,0);
        timeFrame->AddFrame(fTlHSlid);

        fTHSlid = new TGHSlider(timeFrame, 190, kScaleBoth, 1600, kHorizontalFrame, GetDefaultFrameBackground());//,kFALSE, kFALSE, kFALSE, kFALSE);
        fTHSlid->SetRange(0.05, 5.0);
        timeFrame->AddFrame(fTHSlid, new TGLayoutHints(kLHintsExpandX));
        fTHSlid->Associate(this);

        TGHorizontalFrame *fHframe2 = new TGHorizontalFrame(evtidFrame);
        fTeh1 = new TGTextEntry(fHframe2, fTbh1 = new TGTextBuffer(5), 1700);
        fTeh1->SetToolTipText("Time (ns)");
        fTbh1->AddText(0, "0.0");
        fHframe2->AddFrame(fTeh1,new TGLayoutHints(kLHintsExpandX));
        fTeh1->Associate(this);
        TGTextButton *Gobutton         = new TGTextButton(navFrame, "&Go", 1999);
        navFrame->AddFrame(Gobutton, new TGLayoutHints(kLHintsLeft,3,0,3,0));         
        Gobutton->Associate(this);

        //Add Mu2e logo
        std::string logoFile = "TEveEventDisplay/src/Icons/mu2e_logo_oval.png";
        const TGPicture *logo = gClient->GetPicture(logoFile.c_str());
        TGIcon *icon = new TGIcon(navFrame,logo,50,50);
        navFrame->AddFrame(icon,new TGLayoutHints(kLHintsLeft,20,0,0,0));

        celabel = new TGLabel(evtidFrame, "Cluster Energy");

        TGHorizontalFrame *centenergyframe = new TGHorizontalFrame(evtidFrame);
        cminenergy = new TGTextEntry(centenergyframe, _clustminenergy = new TGTextBuffer(5), 1701);
        _clustminenergy->AddText(0, "0.0");
        centenergyframe->AddFrame(cminenergy,new TGLayoutHints(kLHintsExpandX));
        cminenergy->Associate(this);
        spacer = new TGLabel(centenergyframe,"   ");
        centenergyframe->AddFrame(spacer);
        cmaxenergy = new TGTextEntry(centenergyframe, _clustmaxenergy = new TGTextBuffer(5), 1701);
        _clustmaxenergy->AddText(0, "0.0");
        centenergyframe->AddFrame(cmaxenergy,new TGLayoutHints(kLHintsExpandX));
        cmaxenergy->Associate(this);

        std::string clusterenergy = "TEveEventDisplay/src/Icons/purplegradient.png";
        const TGPicture *ce = gClient->GetPicture(clusterenergy.c_str());
        TGIcon *ceicon = new TGIcon(evtidFrame, ce, 40, 8);

        helabel = new TGLabel(evtidFrame, "Hit Energy");
        std::string hitenergy = "TEveEventDisplay/src/Icons/greengradient.png";
        const TGPicture *he = gClient->GetPicture(hitenergy.c_str());
        TGIcon *heicon = new TGIcon(evtidFrame, he, 40 ,8);

        TGHorizontalFrame *hentenergyframe = new TGHorizontalFrame(evtidFrame);
        hminenergy = new TGTextEntry(hentenergyframe, _hitminenergy = new TGTextBuffer(5), 1702);
        _hitminenergy->AddText(0, "0.0");
        hentenergyframe->AddFrame(hminenergy,new TGLayoutHints(kLHintsExpandX));
        hminenergy->Associate(this);
        spacer1 = new TGLabel(hentenergyframe,"   ");
        hentenergyframe->AddFrame(spacer);
        hmaxenergy = new TGTextEntry(hentenergyframe, _hitmaxenergy = new TGTextBuffer(5), 1702);
        _hitmaxenergy->AddText(0, "0.0");
        hentenergyframe->AddFrame(hmaxenergy,new TGLayoutHints(kLHintsExpandX));
        hmaxenergy->Associate(this);

        br = new TGButtonGroup(evtidFrame, "Data Products", kVerticalFrame);
        clusterscheck = new TGCheckButton(br, new TGHotString("Clusters"), 1200);
        clusterscheck->SetState(kButtonDown);
        hitscheck = new TGCheckButton(br, new TGHotString("Hits"), 1201);
        hitscheck->SetState(kButtonDown);
        trackscheck = new TGCheckButton(br, new TGHotString("Tracks"), 1202);
        trackscheck->SetState(kButtonDown);
        cosmicscheck = new TGCheckButton(br, new TGHotString("Cosmics"), 1203);
        cosmicscheck->SetState(kButtonDown);
        cosmictrkscheck = new TGCheckButton(br, new TGHotString("Cosmic Tracks"), 1204);
        cosmictrkscheck->SetState(kButtonDown);
        mctrajcheck = new TGCheckButton(br, new TGHotString("MC Trajectories"), 1205);	
        mctrajcheck->SetState(kButtonDown);
        
        // ... Add horizontal run & event number subframes to vertical evtidFrame
        evtidFrame->AddFrame(runoFrame,new TGLayoutHints(kLHintsExpandX));
        evtidFrame->AddFrame(evnoFrame,new TGLayoutHints(kLHintsExpandX));

        evtidFrame->AddFrame(timeFrame,new TGLayoutHints(kLHintsExpandX));
        evtidFrame->AddFrame(fHframe2,new TGLayoutHints(kLHintsExpandX));
        evtidFrame->AddFrame(celabel, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(ceicon, new TGLayoutHints(kLHintsLeft,20,0,0,0));
        evtidFrame->AddFrame(centenergyframe, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(helabel, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(heicon, new TGLayoutHints(kLHintsLeft,20,0,0,0));
        evtidFrame->AddFrame(hentenergyframe, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(br, new TGLayoutHints(kLHintsExpandX));

        clusterscheck->Associate(this);
        hitscheck->Associate(this);
        trackscheck->Associate(this);
        cosmictrkscheck->Associate(this);
        mctrajcheck->Associate(this);
        
        // ... Add navFrame and evtidFrame to MainFrame
        frmMain->AddFrame(navFrame);
        TGHorizontal3DLine *separator = new TGHorizontal3DLine(frmMain);
        frmMain->AddFrame(separator, new TGLayoutHints(kLHintsExpandX));
        frmMain->AddFrame(evtidFrame);

        frmMain->MapSubwindows();
        frmMain->Resize();
        frmMain->MapWindow();

        browser->StopEmbedding();
        browser->SetTabTitle("Event Nav", 0);
      }
      gEve->AddEvent(new TEveEventManager("Event", "Empty Event"));
      TGLViewer *glv = gEve->GetDefaultGLViewer();
      glv->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
      glv->CurrentCamera().RotateRad(camRotateCenterH_,camRotateCenterV_);
      glv->CurrentCamera().Dolly(camDollyDelta_,kFALSE,kFALSE);
    }

    void TEveMu2eMainWindow::StartCRVProjectionTab(){
      // Create detector and event scenes for ortho views
      CRV2Dproj->fDetXYScene = gEve->SpawnNewScene("CRV Top", "");
      CRV2Dproj->fEvtXYScene = gEve->SpawnNewScene("CRV Top event", "");
    
      // Create XY/RZ calo2Dprojection mgrs, draw projected axes, & add them to scenes
      CRV2Dproj->fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
      TEveProjectionAxes* axes_xy1 = new TEveProjectionAxes(CRV2Dproj->fXYMgr);
      CRV2Dproj->fDetXYScene->AddElement(axes_xy1);
      gEve->AddToListTree(axes_xy1,kTRUE);
      gEve->AddToListTree(CRV2Dproj->fXYMgr,kTRUE);

      // Create side-by-side ortho D1, D2 views in new tab & add det/evt scenes
      TEveWindowSlot *slot = 0;
      TEveWindowPack *pack = 0;

      slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
      pack = slot->MakePack();
      pack->SetElementName("CRV Views");
      pack->SetHorizontal();
      pack->SetShowTitleBar(kFALSE);

      pack->NewSlot()->MakeCurrent();
      CRV2Dproj->fXYView = gEve->SpawnNewViewer("CRV Top View", "");
      CRV2Dproj->fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
      CRV2Dproj->fXYView->AddScene(CRV2Dproj->fDetXYScene);
      CRV2Dproj->fXYView->AddScene(CRV2Dproj->fEvtXYScene);
  }

  void TEveMu2eMainWindow::StartCaloProjectionTab(){
    // Create detector and event scenes for ortho views
    calo2Dproj->fDetXYScene = gEve->SpawnNewScene("Calo XY D0 Scene", "");
    calo2Dproj->fDetRZScene = gEve->SpawnNewScene("Calo XY D1 Scene", "");
    calo2Dproj->fEvtXYScene = gEve->SpawnNewScene("Calo Evt XY D0 Scene", "");
    calo2Dproj->fEvtRZScene = gEve->SpawnNewScene("Calo Evt XY D1 Scene", "");

    // Create XY/RZ calo2Dprojection mgrs, draw projected axes, & add them to scenes
    calo2Dproj->fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    TEveProjectionAxes* axes_xy = new TEveProjectionAxes(calo2Dproj->fXYMgr);
    calo2Dproj->fDetXYScene->AddElement(axes_xy);
    gEve->AddToListTree(axes_xy,kTRUE);
    gEve->AddToListTree(calo2Dproj->fXYMgr,kTRUE);

    calo2Dproj->fRZMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    TEveProjectionAxes* axes_rz = new TEveProjectionAxes(calo2Dproj->fRZMgr);
    calo2Dproj->fDetRZScene->AddElement(axes_rz);
    gEve->AddToListTree(axes_rz,kTRUE);
    gEve->AddToListTree(calo2Dproj->fRZMgr,kTRUE);

    // Create side-by-side ortho D1, D2 views in new tab & add det/evt scenes
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;

    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack();
    pack->SetElementName("Calo Views");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);

    pack->NewSlot()->MakeCurrent();
    calo2Dproj->fXYView = gEve->SpawnNewViewer("Disk0 View", "");
    calo2Dproj->fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    calo2Dproj->fXYView->AddScene(calo2Dproj->fDetXYScene);
    calo2Dproj->fXYView->AddScene(calo2Dproj->fEvtXYScene);

    pack->NewSlot()->MakeCurrent();
    calo2Dproj->fRZView = gEve->SpawnNewViewer("Disk1 View", "");
    calo2Dproj->fRZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    calo2Dproj->fRZView->AddScene(calo2Dproj->fDetRZScene);
    calo2Dproj->fRZView->AddScene(calo2Dproj->fEvtRZScene);

    gEve->GetBrowser()->GetTabRight()->SetTab(0);
  }

  void TEveMu2eMainWindow::StartTrackerProjectionTab(){
    // Create detector and event scenes for ortho views
    tracker2Dproj->fDetXYScene = gEve->SpawnNewScene("Tracker Det XY Scene", "");
    tracker2Dproj->fDetRZScene = gEve->SpawnNewScene("Tracker Det RZ Scene", "");
    tracker2Dproj->fEvtXYScene = gEve->SpawnNewScene("Tracker Evt XY Scene", "");
    tracker2Dproj->fEvtRZScene = gEve->SpawnNewScene("Tracker Evt RZ Scene", "");

    // Create XY/RZ tracker2Dprojection mgrs, draw projected axes, & add them to scenes
    tracker2Dproj->fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    TEveProjectionAxes* axes_xy = new TEveProjectionAxes(tracker2Dproj->fXYMgr);
    tracker2Dproj->fDetXYScene->AddElement(axes_xy);
    gEve->AddToListTree(axes_xy,kTRUE);
    gEve->AddToListTree(tracker2Dproj->fXYMgr,kTRUE);

    tracker2Dproj->fRZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
    TEveProjectionAxes* axes_rz = new TEveProjectionAxes(tracker2Dproj->fRZMgr);
    tracker2Dproj->fDetRZScene->AddElement(axes_rz);
    gEve->AddToListTree(axes_rz,kTRUE);
    gEve->AddToListTree(tracker2Dproj->fRZMgr,kTRUE);

    // Create side-by-side ortho XY & RZ views in new tab & add det/evt scenes
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;

    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack();
    pack->SetElementName("Tracker Views");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);

    pack->NewSlot()->MakeCurrent();
    tracker2Dproj->fXYView = gEve->SpawnNewViewer("Tracker XY View", "");
    tracker2Dproj->fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    tracker2Dproj->fXYView->AddScene(tracker2Dproj->fDetXYScene);
    tracker2Dproj->fXYView->AddScene(tracker2Dproj->fEvtXYScene);

    pack->NewSlot()->MakeCurrent();
    tracker2Dproj->fRZView = gEve->SpawnNewViewer("Tracker RZ View", "");
    tracker2Dproj->fRZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    tracker2Dproj->fRZView->AddScene(tracker2Dproj->fDetRZScene);
    tracker2Dproj->fRZView->AddScene(tracker2Dproj->fEvtRZScene);

    gEve->GetBrowser()->GetTabRight()->SetTab(0);

  }

  void TEveMu2eMainWindow::PrepareCaloProjectionTab(const art::Run& run){
    calo2Dproj->fDetXYScene->DestroyElements();
    calo2Dproj->fDetRZScene->DestroyElements();
    TEveElementList *orthodet0 = new TEveElementList("CaloOrthoDet0");
    TEveElementList *orthodet1 = new TEveElementList("CaloOrthoDet1");
    TGeoVolume* topvol = geom->GetTopVolume(); 
    Mu2eCalo->DrawCaloDetector(run, topvol,orthodet0,orthodet1);

    gEve->AddGlobalElement(orthodet0);
    gEve->AddGlobalElement(orthodet1);

    // ... Import elements of the list into the projected views
    calo2Dproj->fXYMgr->ImportElements(orthodet0, calo2Dproj->fDetXYScene);
    calo2Dproj->fRZMgr->ImportElements(orthodet1, calo2Dproj->fDetRZScene);

    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("OrthoDet")->SetRnrState(kFALSE);

    // ... Turn ON rendering of detector in RPhi and RZ views
    calo2Dproj->fDetXYScene->FindChild("CaloOrthoDet0 [P]")->SetRnrState(kTRUE);
    calo2Dproj->fDetRZScene->FindChild("CaloOrthoDet1 [P]")->SetRnrState(kTRUE);
  }

  void TEveMu2eMainWindow::PrepareTrackerProjectionTab(const art::Run& run){
    tracker2Dproj->fDetXYScene->DestroyElements();
    tracker2Dproj->fDetRZScene->DestroyElements();
    TEveElementList *orthodet = new TEveElementList("OrthoDet");
    TGeoVolume* topvol = geom->GetTopVolume();
    Mu2eTracker->DrawTrackerDetector(run, topvol, orthodet);

    gEve->AddGlobalElement(orthodet);

    // ... Import elements of the list into the projected views
    tracker2Dproj->fXYMgr->ImportElements(orthodet, tracker2Dproj->fDetXYScene);
    tracker2Dproj->fRZMgr->ImportElements(orthodet, tracker2Dproj->fDetRZScene);

    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("OrthoDet")->SetRnrState(kFALSE);

    // ... Turn ON rendering of detector in RPhi and RZ views
    tracker2Dproj->fDetXYScene->FindChild("OrthoDet [P]")->SetRnrState(kTRUE);
    tracker2Dproj->fDetRZScene->FindChild("OrthoDet [P]")->SetRnrState(kTRUE);
  }

  void TEveMu2eMainWindow::PrepareCRVProjectionTab(const art::Run& run){

  CRV2Dproj->fDetXYScene->DestroyElements();
  CRV2Dproj->fEvtXYScene->DestroyElements();

  TGeoVolume* topvol = geom->GetTopVolume();

  TEveElementList *orthodetT1 = new TEveElementList("CRVT1OrthoDet");
  TEveElementList *orthodetT2 = new TEveElementList("CRVT2OrthoDet");
  TEveElementList *orthodetT3 = new TEveElementList("CRVT3OrthoDet");
  TEveElementList *orthodetT4 = new TEveElementList("CRVT4OrthoDet");
  TEveElementList *orthodetlist[] = {orthodetT1, orthodetT2, orthodetT3, orthodetT4};

  Mu2eCRV->DrawCRVDetector(run, topvol, orthodetlist);

  for (unsigned int i=0; i<4; i++){
    gEve->AddGlobalElement(orthodetlist[i]);
  }

  for (unsigned int i=0; i<4; i++){
    CRV2Dproj->fXYMgr->ImportElements(orthodetlist[i], CRV2Dproj->fDetXYScene);
  }

  // ... Turn OFF rendering of duplicate detector in main 3D view
  gEve->GetGlobalScene()->FindChild("OrthoDet")->SetRnrState(kFALSE);

  // ... Turn ON rendering of detector in RPhi and RZ views
  //CRV2Dproj->fDetXYScene->FindChild("OrthoDet [P]")->SetRnrState(kTRUE);
  //CRV2Dproj->fDetRZScene->FindChild("OrthoDets0 [P]")->SetRnrState(kTRUE);

  }

  void TEveMu2eMainWindow::SetRunGeometry(const art::Run& run, int _diagLevel, bool _showBuilding, bool _showDSOnly, bool _showCRV){
    if(gGeoManager){
      gGeoManager->GetListOfNodes()->Delete();
      gGeoManager->GetListOfVolumes()->Delete();
      gGeoManager->GetListOfShapes()->Delete();
    }
    gEve->GetGlobalScene()->DestroyElements();

    // Import the GDML of entire Mu2e Geometry
    geom = mu2e_geom->Geom_Interface::getGeom("TEveEventDisplay/src/fix.gdml");

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

    setRecursiveColorTransp(etopnode->GetNode()->GetVolume(), kWhite-10,70);

    if(!_showBuilding){   
      mu2e_geom->SolenoidsOnly(topnode);
      mu2e_geom->hideTop(topnode, _diagLevel);
    }
    if(_showDSOnly) mu2e_geom->InsideDS(topnode, false );
  
    if(_showCRV) mu2e_geom->InsideCRV(topnode, true);

    //Add static detector geometry to global scene
    gEve->AddGlobalElement(etopnode);
    geom->Draw("ogl");
  }

  void TEveMu2eMainWindow::RedrawDataProducts(std::string type){
    if (type == "Clusters"){
      *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _emptydata.clustercol, calo2Dproj, texttime, true, _show2D, fclustmin, fclustmax, _accumulate);
    }
    if (type == "Hits"){
      if (_data.chcol !=0){*hitenergy = pass_data->AddComboHits(_firstLoop, _emptydata.chcol, tracker2Dproj, texttime, true, _show2D, fhitmin, fhitmax,_accumulate);}
      if(_data.cryHitcol !=0){pass_data->AddCrystalHits(_firstLoop, _emptydata.cryHitcol, calo2Dproj, texttime, true, _show2D,_accumulate);}
    }
    if (type == "Tracks"){
      pass_data->AddHelixPieceWise(_firstLoop, _emptydata.kalseedcol, tracker2Dproj, texttime, true, _show2D, _accumulate);
    }
    if (type == "Cosmics"){
      if(_data.crvcoincol!= 0){pass_data->AddCRVInfo(_firstLoop, _emptydata.crvcoincol, texttime, true, _show2D,  _accumulate);}
      }
    if (type == "Cosmic Tracks"){
      if(_data.cosmiccol!=0){pass_data->AddCosmicTrack(_firstLoop, _emptydata.cosmiccol, tracker2Dproj, texttime, true, _show2D, _accumulate);}
    }
    if (type == "MC Trajectories"){
      if(_data.mctrajcol!=0){pass_mc->AddMCTrajectory(_firstLoop, _emptydata.mctrajcol, tracker2Dproj, true, _show2D, _accumulate);}
    }
    gSystem->ProcessEvents();
    gClient->NeedRedraw(fTeRun);
    gApplication->Run(true);
	}

  void TEveMu2eMainWindow::RedrawGeometry(){
    geom = mu2e_geom->Geom_Interface::getGeom("TEveEventDisplay/src/fix.gdml");

    TGeoNode* topnode = gGeoManager->GetTopNode();
    TEveGeoTopNode* etopnode = new TEveGeoTopNode(gGeoManager, topnode);
    etopnode->SetVisLevel(4);
    etopnode->GetNode()->GetVolume()->SetVisibility(kFALSE);

    setRecursiveColorTransp(etopnode->GetNode()->GetVolume(), kWhite-10,70);
    mu2e_geom->InsideCRV(topnode, true);

    //Add static detector geometry to global scene
    gEve->AddGlobalElement(etopnode);
    geom->Draw("ogl");
  }
  
  Bool_t TEveMu2eMainWindow::ProcessMessage(Long_t msg, Long_t param1, Long_t param2){
    switch (GET_MSG(msg))
    {  
    case kC_HSLIDER:
      if(param1==1600){
        char buf[32]; 
        sprintf(buf, "%.3d", fTHSlid->GetPosition());
        fTbh1->Clear();
        fTbh1->AddText(0, buf);
        fTeh1->SetCursorPosition(fTeh1->GetCursorPosition());
        fTeh1->Deselect();
        gClient->NeedRedraw(fTeh1);
        texttime = fTHSlid->GetPosition();
        pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, texttime, false, _show2D, _accumulate);//, CRV2Dproj);
        *hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, texttime, false, _show2D, fhitmin, fhitmax, _accumulate);
        *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj, texttime, false, _show2D, fclustmin, fclustmax, _accumulate);
        pass_data->AddHelixPieceWise(_firstLoop, _data.kalseedcol,tracker2Dproj, texttime, false, _show2D, _accumulate);
	    }
     break; 
    case kC_TEXTENTRY:
    switch (GET_SUBMSG(msg)){
      case kTE_TEXTCHANGED:
      if (param1 == 1700){
        fTHSlid->SetPosition(atof(fTbh1->GetString()));
      }	
      if (param1 == 1701){
        fclustmin = atof(_clustminenergy->GetString());
        fclustmax = atof(_clustmaxenergy->GetString());
        if (fclustmin < fclustmax) {*clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj, texttime, false, _show2D, fclustmin, fclustmax, _accumulate);}
      }
      if (param1 == 1702){
        fhitmin = atof(_hitminenergy->GetString());
        fhitmax = atof(_hitmaxenergy->GetString());
        if (fhitmin < fhitmax) {*hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, texttime, false, _show2D, fhitmin, fhitmax, _accumulate);}
	    }	
      break;
    }
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
      case kCM_CHECKBUTTON:
      if(param1==1200){
        if(clusterscheck->IsDown()){
        *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj, texttime, false, _show2D, fclustmin, fclustmin, _accumulate);}
        if(!clusterscheck->IsDown() && _data.clustercol!=0){RedrawDataProducts("Clusters");}
      }
      if(param1==1201){
        if(hitscheck->IsDown()){
        *hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, texttime, false, _show2D, fhitmin, fhitmax, _accumulate);
        pass_data->AddCrystalHits(_firstLoop, _data.cryHitcol, calo2Dproj, texttime, false, _show2D, _accumulate);
        }
        if(!hitscheck->IsDown()){RedrawDataProducts("Hits");}
      }
      if(param1==1202){
        if(trackscheck->IsDown()){pass_data->AddHelixPieceWise(_firstLoop, _data.kalseedcol, tracker2Dproj, texttime, false, _show2D, _accumulate);}
        if(!trackscheck->IsDown() && _data.kalseedcol!=0){RedrawDataProducts("Tracks");}
      }
      if(param1==1203){
        if(cosmicscheck->IsDown()){
        pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, texttime, false, _show2D, _accumulate);
        
        }
        if(!cosmicscheck->IsDown()){RedrawDataProducts("Cosmics");}
      }
      if(param1==1204){
        if(cosmictrkscheck->IsDown()){
          pass_data->AddCosmicTrack(_firstLoop, _data.cosmiccol,  tracker2Dproj, texttime, false, _show2D, _accumulate);		
          }
        if(!cosmictrkscheck->IsDown()){RedrawDataProducts("Cosmic Tracks");}
      }
      if(param1==1205){
        if(mctrajcheck->IsDown()){
          pass_mc->AddMCTrajectory(_firstLoop, _data.mctrajcol, tracker2Dproj, false, _show2D, _accumulate);	
        }
        if(!mctrajcheck->IsDown()){RedrawDataProducts("MC Trajectories");}
      }
      break;
      case kCM_BUTTON: 
       if(param1==1111)
       {}
       if(param1==1001)//Forward
       {
         gApplication->Terminate(0);
       }
       if(param1==1100)//Back
       {
          std::cout<<"Still developing backwards navigation"<<std::endl;
       }
       if(param1==1999)//Go
       {
          eventToFind = atoi(fTeEvt->GetText());
          runToFind = atoi(fTeRun->GetText());
          gApplication->Terminate();
       } 
       if(param1==1400){
          RedrawGeometry(); 
      }
      break;
	  } 
  break;
  }
  return kTRUE;
  }

  void TEveMu2eMainWindow::setEvent(const art::Event& event, bool firstLoop, Data_Collections &data, double time, bool show2D, bool accumulate)
  {
    _event=event.id().event();
    _subrun=event.id().subRun();
    _run=event.id().run();
    _firstLoop = firstLoop;
    _show2D = show2D;
    _accumulate = accumulate;
    _data.chcol = data.chcol; 
    _data.clustercol = data.clustercol;
    _data.crvcoincol = data.crvcoincol;
    _data.kalseedcol = data.kalseedcol;
    _data.mctrajcol = data.mctrajcol;
    _data.crvcoincol = data.crvcoincol;
    _data.cryHitcol = data.cryHitcol;
    _data.cosmiccol = data.cosmiccol;
    if (texttime == -1){
      std::vector<double> times = pass_data->getTimeRange(firstLoop, data.chcol, data.crvcoincol, data.clustercol, data.cryHitcol);
      fTHSlid->SetRange(times.at(0), times.at(1));
    }
    if(_data.crvcoincol!=0) pass_data->AddCRVInfo(firstLoop, data.crvcoincol, time, false, show2D, _accumulate);
    hitenergy = new vector<double>(2);
    if(_data.chcol!=0) *hitenergy = pass_data->AddComboHits(firstLoop, data.chcol, tracker2Dproj, time, false, show2D, fhitmin, fhitmax, _accumulate);

    clusterenergy = new vector<double>(2);
    if(_data.clustercol!=0) *clusterenergy = pass_data->AddCaloClusters(firstLoop, data.clustercol, calo2Dproj, time, false, show2D, fclustmin, fclustmax, _accumulate);
    if (_data.cryHitcol!=0) pass_data->AddCrystalHits(_firstLoop, _data.cryHitcol, calo2Dproj, texttime, false, _show2D, _accumulate);
    if(_data.kalseedcol!=0) pass_data->AddHelixPieceWise(firstLoop, data.kalseedcol, tracker2Dproj, time, false, show2D, _accumulate);
    if(_data.cosmiccol!=0) pass_data->AddCosmicTrack(firstLoop, data.cosmiccol, tracker2Dproj, time, false, show2D, _accumulate);
    if(_data.mctrajcol!=0) pass_mc->AddMCTrajectory(firstLoop, data.mctrajcol, tracker2Dproj, false, show2D, _accumulate);
    
    gSystem->ProcessEvents();
    gSystem->IgnoreInterrupt();
    gSystem->IgnoreSignal(kSigTermination);
    gSystem->IgnoreSignal(kSigSegmentationViolation);
 
    gClient->NeedRedraw(fTeRun);
    _clustminenergy->Clear();
    _clustmaxenergy->Clear();
    _hitminenergy->Clear();
    _hitmaxenergy->Clear();
    _clustminenergy->AddText(0, (to_string(clusterenergy->at(0))).c_str());
    _clustmaxenergy->AddText(0, (to_string(clusterenergy->at(1))).c_str());

    _hitminenergy->AddText(0, (to_string(hitenergy->at(0))).c_str());
    _hitmaxenergy->AddText(0, (to_string(hitenergy->at(1))).c_str());
    gApplication->Run(true);

    gEve->Redraw3D(kTRUE);
  }

  void TEveMu2eMainWindow::fillEvent(bool firstLoop)
   {
    std::string eventInfoText;
    eventInfoText=Form("Event #: %i",_event);
    if(_eventNumberText==nullptr) 
    {
      _eventNumberText = new TText(0.6,-0.8, eventInfoText.c_str());
      _eventNumberText->SetTextColor(5);
      _eventNumberText->SetTextSize(0.025);
      _eventNumberText->Draw("same");
    }
    else _eventNumberText->SetTitle(eventInfoText.c_str());
    eventInfoText=Form("Sub Run #: %i",_subrun);
    if(_subrunNumberText==nullptr)
    {
      _subrunNumberText = new TText(0.6,-0.75,eventInfoText.c_str());
      _subrunNumberText->SetTextColor(5);
      _subrunNumberText->SetTextSize(0.025);
      _subrunNumberText->Draw("same");
    }
    else _subrunNumberText->SetTitle(eventInfoText.c_str());
    eventInfoText=Form("Run #: %i",_run);
    if(_runNumberText==nullptr)
    {
      _runNumberText = new TText(0.6,-0.7,eventInfoText.c_str());
      _runNumberText->SetTextColor(5);
      _runNumberText->SetTextSize(0.025);
      _runNumberText->Draw("same");
    }
    else _runNumberText->SetTitle(eventInfoText.c_str());
    if(_timeText==nullptr)
    {
      _timeText = new TText(0.6,-0.7,eventInfoText.c_str());
      _timeText->SetTextColor(5);
      _timeText->SetTextSize(0.025);
      _timeText->Draw("same");
    }
    else _timeText->SetTitle(eventInfoText.c_str());
    this->Layout();
  }

  int TEveMu2eMainWindow::getEventToFind(bool &findEvent) const
  {
    findEvent=_findEvent;
    return _eventToFind;
  }

  bool TEveMu2eMainWindow::isClosed() const
  {
    return _isClosed;
  }
}
