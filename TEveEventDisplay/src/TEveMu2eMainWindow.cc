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
#include<TGLEmbeddedViewer.h>
#include <TGMsgBox.h>
#include<TPolyLine3D.h>
#include <TGSplitFrame.h>
// ... libRGL
#include <TGLViewer.h>
#include <TVirtualX.h>
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
    //fTHSlid(0),
    fTlRun(0),
    fTlEvt(0),
    fTlTEvt(0),
    //fTlHSlid(0),
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
      browser = gEve->GetBrowser();
      CreateGUI();
      CreateMultiViews();
      gEve->AddEvent(new TEveEventManager("Event", "Empty Event"));

    }

  void TEveMu2eMainWindow::CreateMultiViews(){
   gEve->GetBrowser()->GetTabRight()->SetTab(0);

   fPad = new TEvePad();
   fPad->SetFillColor(kBlack);
   // create the split frames
   fSplitFrame = new TGSplitFrame(this, 900, 1300);
   AddFrame(fSplitFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   // split it once
   fSplitFrame->HSplit(350);
   fSplitFrame->GetFirst()->VSplit(410);
   fSplitFrame->GetSecond()->VSplit(410);
  // get top (main) split frame


   frm = fSplitFrame->GetFirst()->GetFirst();
   frm->SetName("Calorimeter_XY_View");
   fViewer0 = new TGLEmbeddedViewer(frm, fPad);
   frm->AddFrame(fViewer0->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer0->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[0] = new TEveViewer("SplitGLViewer[0]");
   fViewer[0]->SetGLViewer(fViewer0, fViewer0->GetFrame());
   fViewer[0]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[0]);
     proj0 = gEve->SpawnNewScene("Calorimeter XY Scene");
     //fViewer[1]->AddScene(fdetXY);
     CfXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
     proj0->AddElement(CfXYMgr);
     TEveProjectionAxes* axes_xy = new TEveProjectionAxes(CfXYMgr);
     proj0->AddElement(axes_xy);
     gEve->AddToListTree(axes_xy,kTRUE);
     gEve->AddToListTree(CfXYMgr,kTRUE);
     fViewer[0]->AddScene(proj0);
  }

   frm = fSplitFrame->GetFirst()->GetSecond();
   frm->SetName("Calorimeter_RZ_View");
   fViewer1 = new TGLEmbeddedViewer(frm, fPad);
   frm->AddFrame(fViewer1->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer1->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[1] = new TEveViewer("SplitGLViewer[0]");
   fViewer[1]->SetGLViewer(fViewer1, fViewer1->GetFrame());
   fViewer[1]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[1]);
     proj1 = gEve->SpawnNewScene("Calorimeter XY Scene");
     //fViewer[1]->AddScene(fdetXY);
     CfRZMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
     proj1->AddElement(CfRZMgr);
     TEveProjectionAxes* axes_xy = new TEveProjectionAxes(CfRZMgr);
     proj1->AddElement(axes_xy);
     gEve->AddToListTree(axes_xy,kTRUE);
     gEve->AddToListTree(CfRZMgr,kTRUE);
     fViewer[1]->AddScene(proj1);
}

   frm = fSplitFrame->GetSecond()->GetFirst();
   frm->SetName("Tracker_XY_View");
   fViewer2 = new TGLEmbeddedViewer(frm, fPad);
   frm->AddFrame(fViewer2->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer2->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[2] = new TEveViewer("SplitGLViewer[2]");
   fViewer[2]->SetGLViewer(fViewer2, fViewer2->GetFrame());
   fViewer[2]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[2]);
     proj2 = gEve->SpawnNewScene("Tracker XY Scene");
     //fViewer[1]->AddScene(fdetXY);
     TfXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
     proj2->AddElement(TfXYMgr);
     TEveProjectionAxes* axes_xytracker = new TEveProjectionAxes(TfXYMgr);
     proj2->AddElement(axes_xytracker);
     gEve->AddToListTree(axes_xytracker,kTRUE);
     gEve->AddToListTree(TfXYMgr,kTRUE);
     fViewer[2]->AddScene(proj2);
}

   frm = fSplitFrame->GetSecond()->GetSecond();
   frm->SetName("Tracker_RZ_View");
   fViewer3 = new TGLEmbeddedViewer(frm, fPad);
   frm->AddFrame(fViewer3->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer3->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[3] = new TEveViewer("SplitGLViewer[3]");
   fViewer[3]->SetGLViewer(fViewer3, fViewer3->GetFrame());
   fViewer[3]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[3]);
     proj3 = gEve->SpawnNewScene("Tracker XY Scene");
     //fViewer[1]->AddScene(fdetXY);
     TfRZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
     proj3->AddElement(TfRZMgr);
     TEveProjectionAxes* axes_xytracker = new TEveProjectionAxes(TfRZMgr);
     proj3->AddElement(axes_xytracker);
     gEve->AddToListTree(axes_xytracker,kTRUE);
     gEve->AddToListTree(TfRZMgr,kTRUE);
     fViewer[3]->AddScene(proj3);
}
	
   Resize(GetDefaultSize());
   MapSubwindows();
   MapWindow();

}

  void TEveMu2eMainWindow::CreateGUI(){
      FontStruct_t buttonfont = gClient->GetFontByName("-*-helvetica-medium-r-*-*-8-*-*-*-*-*-iso8859-1");
      GCValues_t gval;
      gval.fMask = kGCForeground | kGCFont;
      gval.fFont = gVirtualX->GetFontHandle(buttonfont);
      gClient->GetColorByName("black", gval.fForeground);

      browser->StartEmbedding(TRootBrowser::kLeft); // insert nav frame as new tab in left pane

      frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
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

        // ... Create run num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* runoFrame = new TGHorizontalFrame(evtidFrame);
        fTlRun = new TGLabel(runoFrame,"Run Number");
        fTlRun->SetTextJustify(kTextLeft);
        fTlRun->SetMargins(5,5,5,0);
        runoFrame->AddFrame(fTlRun);

        fTeRun = new TGTextEntry(runoFrame, _runNumber = new TGTextBuffer(5), 1);
        _runNumber->AddText(0, "Enter Run Number");

        runoFrame->AddFrame(fTeRun,new TGLayoutHints(kLHintsExpandX));
        fTeRun->Associate(this);

        // ... Create evt num text entry widget and connect to "GotoEvent" rcvr in visutils
        TGHorizontalFrame* evnoFrame = new TGHorizontalFrame(evtidFrame);
        fTlEvt = new TGLabel(evnoFrame,"Evt Number");
        fTlEvt->SetTextJustify(kTextLeft);
        fTlEvt->SetMargins(5,5,5,0);
        evnoFrame->AddFrame(fTlEvt);

        fTeEvt = new TGTextEntry(evnoFrame, _eventNumber = new TGTextBuffer(5), 1);
        _eventNumber->AddText(0, "Enter Run Number");

        evnoFrame->AddFrame(fTeEvt,new TGLayoutHints(kLHintsExpandX));
        fTeEvt->Associate(this);
        
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
        timelabel = new TGLabel(evtidFrame, "Time Interval [ns]");
        TGHorizontalFrame *henttimeframe = new TGHorizontalFrame(evtidFrame);
        hmintime = new TGTextEntry(henttimeframe, _hitmintime = new TGTextBuffer(5), 1703);
        _hitmintime->AddText(0, "0.0");
        henttimeframe->AddFrame(hmintime,new TGLayoutHints(kLHintsExpandX));
        hmintime->Associate(this);
        spacer1 = new TGLabel(henttimeframe,"   ");
        henttimeframe->AddFrame(spacer);
        hmaxtime = new TGTextEntry(henttimeframe, _hitmaxtime = new TGTextBuffer(5), 1703);
        _hitmaxtime->AddText(0, "0.0");
        henttimeframe->AddFrame(hmaxtime,new TGLayoutHints(kLHintsExpandX));
        hmaxtime->Associate(this);

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

        evtidFrame->AddFrame(celabel, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(ceicon, new TGLayoutHints(kLHintsLeft,20,0,0,0));
        evtidFrame->AddFrame(centenergyframe, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(helabel, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(heicon, new TGLayoutHints(kLHintsLeft,20,0,0,0));
        evtidFrame->AddFrame(hentenergyframe, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(timelabel, new TGLayoutHints(kLHintsLeft,3,0,3,0));
        evtidFrame->AddFrame(henttimeframe, new TGLayoutHints(kLHintsLeft,3,0,3,0));
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


}

   void TEveMu2eMainWindow::StartProjectionTabs(){
	//pass_proj->CreateCRVProjection(CRV2Dproj);
	pass_proj->CreateCaloProjection(calo2Dproj);
	pass_proj->CreateTrackerProjection(tracker2Dproj);
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


    CfXYMgr->ImportElements(orthodet0);
    CfRZMgr->ImportElements(orthodet1);
    // ... Import elements of the list into the projected views
    calo2Dproj->fXYMgr->ImportElements(orthodet0, calo2Dproj->fDetXYScene);
    calo2Dproj->fRZMgr->ImportElements(orthodet1, calo2Dproj->fDetRZScene);

    //fXYMgr->ImportElements(orthodet0, fdetXY);
    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("OrthoDet")->SetRnrState(kFALSE);

    // ... Turn ON rendering of detector in RPhi and RZ views
    calo2Dproj->fDetXYScene->FindChild("CaloOrthoDet0 [P]")->SetRnrState(kTRUE);
    calo2Dproj->fDetRZScene->FindChild("CaloOrthoDet1 [P]")->SetRnrState(kTRUE);


  }

  void TEveMu2eMainWindow::PrepareTrackerProjectionTab(const art::Run& run){
    tracker2Dproj->fDetXYScene->DestroyElements();
    tracker2Dproj->fDetRZScene->DestroyElements();
    //fdetXY->DestroyElements();

    TEveElementList *orthodet = new TEveElementList("OrthoDet");
    TEveElementList *orthodetsplit = new TEveElementList("OrthoDet");
    TGeoVolume* topvol = geom->GetTopVolume();
    Mu2eTracker->DrawTrackerDetector(run, topvol, orthodet);
    Mu2eTracker->DrawTrackerDetector(run, topvol, orthodetsplit);

    gEve->AddGlobalElement(orthodet);


    // ... Import elements of the list into the projected views

    TfXYMgr->ImportElements(orthodetsplit);
    TfRZMgr->ImportElements(orthodetsplit);

    tracker2Dproj->fXYMgr->ImportElements(orthodet, tracker2Dproj->fDetXYScene);	
    tracker2Dproj->fRZMgr->ImportElements(orthodet, tracker2Dproj->fDetRZScene);


    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("OrthoDet")->SetRnrState(kFALSE);


    //fdetXY->FindChild("OrthoDet [P]")->SetRnrState(kTRUE);
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
      *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _emptydata.clustercol, calo2Dproj,true, fclustmin, fclustmax, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
    }
    if (type == "Hits"){
      if (_data.chcol !=0){*hitenergy = pass_data->AddComboHits(_firstLoop, _emptydata.chcol, tracker2Dproj,  true, fhitmin, fhitmax,ftimemin, ftimemax,_accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
      if(_data.cryHitcol !=0){pass_data->AddCrystalHits(_firstLoop, _emptydata.cryHitcol, calo2Dproj, ftimemin, ftimemax, true, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);}
    }
    if (type == "Tracks"){
      pass_data->AddHelixPieceWise3D(_firstLoop, _emptydata.track_tuple, tracker2Dproj, ftimemin, ftimemax, true,  _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);

    }
    if (type == "Cosmics"){
      if(_data.crvcoincol!= 0){pass_data->AddCRVInfo(_firstLoop, _emptydata.crvcoincol, ftimemin, ftimemax, true,  _accumulate);}
      }
    if (type == "Cosmic Tracks"){
      if(_data.cosmiccol!=0){pass_data->AddCosmicTrack(_firstLoop, _emptydata.cosmiccol, tracker2Dproj, ftimemin, ftimemax, true, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
    }
    if (type == "MC Trajectories"){
      if(_data.mctrajcol!=0){pass_mc->AddFullMCTrajectory(_firstLoop, _emptydata.mctrajcol, tracker2Dproj, true, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
    }
    gSystem->ProcessEvents();
    gClient->NeedRedraw(fTeRun);
    gApplication->Run(true);
	}

  
  Bool_t TEveMu2eMainWindow::ProcessMessage(Long_t msg, Long_t param1, Long_t param2){
    switch (GET_MSG(msg))
    {  
    
    case kC_TEXTENTRY:
    switch (GET_SUBMSG(msg)){
      case kTE_TEXTCHANGED:
      if (param1 == 1700){
        fTHSlid->SetPosition(atof(fTbh1->GetString()));
      }	
      if (param1 == 1701){
        fclustmin = atof(_clustminenergy->GetString());
        fclustmax = atof(_clustmaxenergy->GetString());
        if (fclustmin < fclustmax) {*clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj,  false, fclustmin, fclustmax, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);}
        if (fclustmin > fclustmax){
          std::cout<<"Cluster Minimum Energy is greater than Maximum Energy"<<std::endl;
          char msg[300];
          sprintf(msg, "Error #%i : Cluster minimum energy larger than maximum", true);
          new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), "Event Not Found", msg, kMBIconExclamation,kMBOk);
        }
      }
      if (param1 == 1702){
        fhitmin = atof(_hitminenergy->GetString());
        fhitmax = atof(_hitmaxenergy->GetString());
        if (fhitmin < fhitmax) {*hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
        if (fhitmin > fhitmax){
          std::cout<<"Hit Minimum Energy is greater than Maximum Energy"<<std::endl;
          char msg[300];
          sprintf(msg, "Error #%i : Hit minimum energy larger than maximum", true);
          new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(), "Event Not Found", msg, kMBIconExclamation,kMBOk);
        }
	    }	
	    if (param1 == 1703){
        ftimemin = atof(_hitmintime->GetString());
        ftimemax = atof(_hitmaxtime->GetString());
        if (ftimemin < ftimemax) {
           if(_data.chcol!=0) {
            *hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, false, fhitmin, fhitmax, ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
           }  if(_data.clustercol!=0) {
            *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj,  false, fclustmin, fclustmax,ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
            }if(_data.crvcoincol!=0) {
               pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, ftimemin, ftimemax, false, _accumulate);
            }
          }
	    }
      break;
    }
    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
      case kCM_CHECKBUTTON:
      if(param1==1200){
        if(clusterscheck->IsDown()){
        *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj, false, fclustmin, fclustmin, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);}
        if(!clusterscheck->IsDown() && _data.clustercol!=0){RedrawDataProducts("Clusters");}
      }
      if(param1==1201){
        if(hitscheck->IsDown()){
        *hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
        pass_data->AddCrystalHits(_firstLoop, _data.cryHitcol, calo2Dproj, ftimemin, ftimemax, false, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
        }
        if(!hitscheck->IsDown()){RedrawDataProducts("Hits");}
      }
      if(param1==1202){

        if(trackscheck->IsDown()){pass_data->AddHelixPieceWise3D(_firstLoop, _data.track_tuple, tracker2Dproj, ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
       
        if(!trackscheck->IsDown()){RedrawDataProducts("Tracks");}
      }
      if(param1==1203){
        if(cosmicscheck->IsDown()){
        pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, ftimemin, ftimemax, false, _accumulate);
        
        }
        if(!cosmicscheck->IsDown()){RedrawDataProducts("Cosmics");}
      }
      if(param1==1204){
        if(cosmictrkscheck->IsDown()){
          pass_data->AddCosmicTrack(_firstLoop, _data.cosmiccol,  tracker2Dproj,  ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);		
          }
        if(!cosmictrkscheck->IsDown()){RedrawDataProducts("Cosmic Tracks");}
      }
      if(param1==1205){
        if(mctrajcheck->IsDown()){
          pass_mc->AddFullMCTrajectory(_firstLoop, _data.mctrajcol, tracker2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);	
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
          usereventSelected = true;
          gApplication->Terminate(0);
       } 
       if(param1==1400){
          //RedrawGeometry(); 
      }
      break;
	  } 
  break;
  }
  return kTRUE;
  }

  void TEveMu2eMainWindow::setEvent(const art::Event& event, bool firstLoop, Data_Collections &data, double time, bool accumulate, int &runn, int &eventn, bool &eventSelected)
  {
    _event=event.id().event();
    _subrun=event.id().subRun();
    _run=event.id().run();
    _firstLoop = firstLoop;
    _accumulate = accumulate;
    _data.chcol = data.chcol; 
    _data.clustercol = data.clustercol;
    _data.crvcoincol = data.crvcoincol;
    _data.track_tuple = data.track_tuple;
    _data.mctrajcol = data.mctrajcol;
    _data.crvcoincol = data.crvcoincol;
    _data.cryHitcol = data.cryHitcol;
    _data.cosmiccol = data.cosmiccol;
    std::vector<const KalSeedCollection*> track_list = std::get<1>(data.track_tuple);
    std::cout<<"in main "<<track_list.size()<<std::endl;
    std::vector<double> times = pass_data->getTimeRange(firstLoop, data.chcol, data.crvcoincol, data.clustercol, data.cryHitcol);

    if(_data.crvcoincol!=0) pass_data->AddCRVInfo(firstLoop, data.crvcoincol, ftimemin, ftimemax, false, _accumulate);
    hitenergy = new vector<double>(2);
    if(_data.chcol!=0) *hitenergy = pass_data->AddComboHits(firstLoop, data.chcol, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);

    clusterenergy = new std::vector<double>(2);

    if(_data.clustercol!=0) *clusterenergy = pass_data->AddCaloClusters(firstLoop, data.clustercol, calo2Dproj, false, fclustmin, fclustmax, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
    if (_data.cryHitcol!=0) pass_data->AddCrystalHits(firstLoop, data.cryHitcol, calo2Dproj, ftimemin, ftimemax, false, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
    pass_data->AddHelixPieceWise3D(firstLoop, data.track_tuple, tracker2Dproj,  ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
    if(_data.cosmiccol!=0) pass_data->AddCosmicTrack(firstLoop, data.cosmiccol, tracker2Dproj, ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
    if(_data.mctrajcol!=0) pass_mc->AddFullMCTrajectory(firstLoop, data.mctrajcol, tracker2Dproj, false, _accumulate,  TfXYMgr, TfRZMgr, proj2, proj3);

    
    gSystem->ProcessEvents();
    gSystem->IgnoreInterrupt();
    gSystem->IgnoreSignal(kSigTermination);
    gSystem->IgnoreSignal(kSigSegmentationViolation);
 
    gClient->NeedRedraw(fTeRun);
    _clustminenergy->Clear();
    _clustmaxenergy->Clear();
    _hitminenergy->Clear();
    _hitmaxenergy->Clear();
    _hitmintime->Clear();
    _hitmaxtime->Clear();
    _clustminenergy->AddText(0, (to_string(clusterenergy->at(0))).c_str());
    _clustmaxenergy->AddText(0, (to_string(clusterenergy->at(1))).c_str());

    _hitminenergy->AddText(0, (to_string(hitenergy->at(0))).c_str());
    _hitmaxenergy->AddText(0, (to_string(hitenergy->at(1))).c_str());
    _hitmintime->AddText(0, (to_string(times.at(0))).c_str());
    _hitmaxtime->AddText(0, (to_string(times.at(1))).c_str());
    gApplication->Run(true);

    gEve->Redraw3D(kTRUE);
    if(usereventSelected == true){
      eventn = eventToFind;
      runn = runToFind;
      eventSelected = true;
    }
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
