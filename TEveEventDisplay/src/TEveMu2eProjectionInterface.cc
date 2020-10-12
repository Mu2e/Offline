#include "TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eProjectionInterface.h"

using namespace mu2e;
namespace mu2e{
   void TEveMu2eProjectionInterface::CreateCRVProjection(TEveMu2e2DProjection *CRV2Dproj){
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



  void TEveMu2eProjectionInterface::CreateCaloProjection(TEveMu2e2DProjection *calo2Dproj){
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

    
   // gEve->GetBrowser()->GetTabRight()->SetTab(0);
  }

  void TEveMu2eProjectionInterface::CreateTrackerProjection(TEveMu2e2DProjection *tracker2Dproj){
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

    //gEve->GetBrowser()->GetTabRight()->SetTab(0);

  }

/*
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

  }*/
}
