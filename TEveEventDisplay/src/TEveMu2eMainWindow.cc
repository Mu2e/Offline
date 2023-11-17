//TEveMu2e:
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMainWindow.h"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
namespace fhicl
{
  class ParameterSet;
}

using namespace mu2e;
using namespace std;

  /*------------Function to make colour scheme:-------------*/
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


  /*------ Default Constructor ------ */
  TEveMu2eMainWindow::TEveMu2eMainWindow() : TGMainFrame(gClient->GetRoot(), 320, 320){}

  /*------------Function to construct main frame, add buttons and GUI:-------------*/
  TEveMu2eMainWindow::TEveMu2eMainWindow(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet _pset, const DrawOptions drawOpts, const GeomOptions geomOpts) :
    TGMainFrame(p, w, h),
    DrawOpts(drawOpts),
    GeomOpts(geomOpts),
    fTeRun(0),
    fTeEvt(0),
    fTTEvt(0),
    fTlRun(0),
    fTlEvt(0),
    fTlTEvt(0),
    br(0),
    clusterscheck(0),
    hitscheck(0),
    trackscheck(0),
    cosmicscheck(0),
    cosmictrkscheck(0),
    mctrajcheck(0)
    {
      //Create TEve Manager:
      TEveManager::Create();
      //Create Browser :
      gEve->GetBrowser()->GetTabRight()->SetTab(0);
      gClient->GetRoot();
      browser = gEve->GetBrowser();
            //Build GUI (function below)
      CreateGUI();
      //Build Multiple View Window:
      //CreateMultiViews(); --> option for pop up window no deprecated
      //Add your Event:
      gEve->AddEvent(new TEveEventManager("Event", "Empty Event"));

    }

  /*------------Function to create multiple windows for 2D projections:-------------*/
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
     proj3 = gEve->SpawnNewScene("Tracker RZ Scene");
     //fViewer[1]->AddScene(fdetXY);
     TfRZMgr = new TEveProjectionManager(TEveProjection::kPT_ZY);
     proj3->AddElement(TfRZMgr);
     TEveProjectionAxes* axes_xytracker = new TEveProjectionAxes(TfRZMgr);
     proj3->AddElement(axes_xytracker);
     gEve->AddToListTree(axes_xytracker,kTRUE);
     gEve->AddToListTree(TfRZMgr,kTRUE);
     fViewer[3]->AddScene(proj3);
   }
          // create the split frames
   fPadCRV = new TEvePad();
   fPadCRV->SetFillColor(kBlack);

   fSplitFrameCRV = new TGSplitFrame(this, 900, 1300);
   AddFrame(fSplitFrameCRV, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   // split it once
   fSplitFrameCRV->HSplit(350);
   //fSplitFrameCRV->GetFirst()->VSplit(410);
   //fSplitFrameCRV->GetSecond()->VSplit(410);
   //get top (main) split frame
   frmCRV = fSplitFrameCRV->GetFirst();//->GetFirst();
   frmCRV->SetName("CRV_XY_View");
   fViewer4 = new TGLEmbeddedViewer(frmCRV, fPadCRV);
   frmCRV->AddFrame(fViewer4->GetFrame(), new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer4->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[4] = new TEveViewer("SplitGLViewer[4]");
   fViewer[4]->SetGLViewer(fViewer4, fViewer4->GetFrame());
   fViewer[4]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[4]);
     proj4 = gEve->SpawnNewScene("CRV XY Scene");
     //fViewer[1]->AddScene(fdetXY);
     CrfXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
     proj4->AddElement(CrfXYMgr);
     TEveProjectionAxes* axes_xy = new TEveProjectionAxes(CrfXYMgr);
     proj4->AddElement(axes_xy);
     gEve->AddToListTree(axes_xy,kTRUE);
     gEve->AddToListTree(CrfXYMgr,kTRUE);
     fViewer[4]->AddScene(proj4);
   }

   frmCRV = fSplitFrameCRV->GetSecond();//->GetSecond();
   frmCRV->SetName("CRV_YZ_View");
   fViewer5 = new TGLEmbeddedViewer(frmCRV, fPadCRV);
   frmCRV->AddFrame(fViewer5->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer5->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to

   fViewer[5] = new TEveViewer("SplitGLViewer[5]");
   fViewer[5]->SetGLViewer(fViewer5, fViewer5->GetFrame());
   fViewer[5]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
     gEve->GetViewers()->AddElement(fViewer[5]);
     proj5 = gEve->SpawnNewScene("CRV YZ Scene");
     //fViewer[1]->AddScene(fdetXY);
     CrfRZMgr = new TEveProjectionManager(TEveProjection::kPT_ZY);
     proj5->AddElement(CrfRZMgr);
     TEveProjectionAxes* axes_xytracker = new TEveProjectionAxes(CrfRZMgr);
     proj5->AddElement(axes_xytracker);
     gEve->AddToListTree(axes_xytracker,kTRUE);
     gEve->AddToListTree(CrfRZMgr,kTRUE);
     fViewer[5]->AddScene(proj5);
   }

   Resize(GetDefaultSize());
   MapSubwindows();
   MapWindow();

}

  /*------------Function to create GUI buttons etc.:-------------*/
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
        ConfigFileLookupPolicy configFile;
        std::string logoFile =  configFile("Offline/TEveEventDisplay/src/Icons/mu2e_logo_oval.png");
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

        std::string clusterenergy =  configFile("Offline/TEveEventDisplay/src/Icons/purplegradient.png");
        const TGPicture *ce = gClient->GetPicture(clusterenergy.c_str());
        TGIcon *ceicon = new TGIcon(evtidFrame, ce, 40, 8);

        helabel = new TGLabel(evtidFrame, "Hit Energy");
        std::string hitenergy = configFile("Offline/TEveEventDisplay/src/Icons/greengradient.png");;
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

  /*------------Function to create 2D Tabs:-------------*/
  void TEveMu2eMainWindow::StartProjectionTabs(){
          if(DrawOpts.addCRVInfo) {pass_proj->CreateCRVProjection(CRV2Dproj);}
          if(DrawOpts.addClusters or DrawOpts.addCryHits)pass_proj->CreateCaloProjection(calo2Dproj);
          pass_proj->CreateTrackerProjection(tracker2Dproj);
  }

    /*------------Function to add calo 2D projection to display:-------------*/
  void TEveMu2eMainWindow::CreateCaloProjection(){
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

  /*------------Function to add tracker 2D projection to display:-------------*/
  void TEveMu2eMainWindow::CreateTrackerProjection(){
    // Create detector and event scenes for ortho views
    tracker2Dproj->fDetXYScene = gEve->SpawnNewScene("Tracker Det XY Scene", "");
    tracker2Dproj->fDetRZScene = gEve->SpawnNewScene("Tracker Det RZ Scene", "");
    tracker2Dproj->fEvtXYScene = gEve->SpawnNewScene("Tracker Evt XY Scene", "");
    tracker2Dproj->fEvtRZScene = gEve->SpawnNewScene("Tracker Evt RZ Scene", "");

    // Create XY/RZ tracker2Dprojection mgrs, draw projected axes, & add them to scenes
    tracker2Dproj->fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    TEveProjectionAxes* axes_xy = new TEveProjectionAxes(tracker2Dproj->fXYMgr);
    tracker2Dproj->fDetXYScene->AddElement(axes_xy);
    tracker2Dproj->fEvtXYScene->AddElement(axes_xy);
    gEve->AddToListTree(axes_xy,kTRUE);
    gEve->AddToListTree(tracker2Dproj->fXYMgr,kTRUE);

    tracker2Dproj->fRZMgr = new TEveProjectionManager(TEveProjection::kPT_ZY);
    TEveProjectionAxes* axes_rz = new TEveProjectionAxes(tracker2Dproj->fRZMgr);
    tracker2Dproj->fDetRZScene->AddElement(axes_rz);
    tracker2Dproj->fEvtRZScene->AddElement(axes_rz);
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
         /*------------Function to add CRV 2D projection to display:-------------*/
  void TEveMu2eMainWindow::CreateCRVProjection(){
    // Create detector and event scenes for ortho views
    CRV2Dproj->fDetXYScene = gEve->SpawnNewScene("CRV Det XY Scene", "");
    CRV2Dproj->fDetRZScene = gEve->SpawnNewScene("CRV Det RZ Scene", "");
    CRV2Dproj->fEvtXYScene = gEve->SpawnNewScene("CRV Evt XY Scene", "");
    CRV2Dproj->fEvtRZScene = gEve->SpawnNewScene("CRV Evt RZ Scene", "");

    // Create XY/RZ CRV2Dprojection mgrs, draw projected axes, & add them to scenes
    CRV2Dproj->fXYMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
    TEveProjectionAxes* axes_xy = new TEveProjectionAxes(CRV2Dproj->fXYMgr);
    CRV2Dproj->fDetXYScene->AddElement(axes_xy);
    CRV2Dproj->fEvtXYScene->AddElement(axes_xy);
    gEve->AddToListTree(axes_xy,kTRUE);
    gEve->AddToListTree(CRV2Dproj->fXYMgr,kTRUE);

    // Create side-by-side ortho XY & RZ views in new tab & add det/evt scenes
    TEveWindowSlot *slot = 0;
    TEveWindowPack *pack = 0;

    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    pack = slot->MakePack();
    pack->SetElementName("CRV XY View");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);

    pack->NewSlot()->MakeCurrent();
    CRV2Dproj->fXYView = gEve->SpawnNewViewer("CRV XY View", "");
    CRV2Dproj->fXYView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    CRV2Dproj->fXYView->AddScene(CRV2Dproj->fDetXYScene);
    CRV2Dproj->fXYView->AddScene(CRV2Dproj->fEvtXYScene);

    gEve->GetBrowser()->GetTabRight()->SetTab(0);

    CRV2Dproj->fRZMgr = new TEveProjectionManager(TEveProjection::kPT_ZY);
    TEveProjectionAxes* axes_rz = new TEveProjectionAxes(CRV2Dproj->fRZMgr);
    CRV2Dproj->fDetRZScene->AddElement(axes_rz);
    CRV2Dproj->fEvtRZScene->AddElement(axes_rz);
    gEve->AddToListTree(axes_rz,kTRUE);
    gEve->AddToListTree(CRV2Dproj->fRZMgr,kTRUE);

    TEveWindowSlot *slotnew = 0;
    TEveWindowPack *packnew = 0;
    slotnew = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    packnew = slotnew->MakePack();
    packnew->SetElementName("CRV YZ View");
    packnew->SetHorizontal();
    packnew->SetShowTitleBar(kFALSE);
    packnew->NewSlot()->MakeCurrent();
    CRV2Dproj->fRZView = gEve->SpawnNewViewer("CRV YZ View", "");
    CRV2Dproj->fRZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    CRV2Dproj->fRZView->AddScene(CRV2Dproj->fDetRZScene);
    CRV2Dproj->fRZView->AddScene(CRV2Dproj->fEvtRZScene);

    gEve->GetBrowser()->GetTabRight()->SetTab(0);

  }


  /*------------Function to create Calo 2D tab:-------------*/
  void TEveMu2eMainWindow::PrepareCaloProjectionTab(const art::Run& run){
    calo2Dproj->fDetXYScene->DestroyElements();
    calo2Dproj->fDetRZScene->DestroyElements();

    TEveElementList *orthodet0 = new TEveElementList("CaloOrthoDet0");
    TEveElementList *orthodet1 = new TEveElementList("CaloOrthoDet1");
    TGeoVolume* topvol = geom->GetTopVolume();
    Mu2eCalo->DrawCaloDetector(run, topvol,orthodet0,orthodet1);
    gEve->AddGlobalElement(orthodet0);
    gEve->AddGlobalElement(orthodet1);

    //CfXYMgr->ImportElements(orthodet0);
    //CfRZMgr->ImportElements(orthodet1);

    // ... Import elements of the list into the projected views
    calo2Dproj->fXYMgr->ImportElements(orthodet0, calo2Dproj->fDetXYScene);
    calo2Dproj->fRZMgr->ImportElements(orthodet1, calo2Dproj->fDetRZScene);

    //fXYMgr->ImportElements(orthodet0, fdetXY);
    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("CaloOrthoDet0")->SetRnrState(kFALSE);
    gEve->GetGlobalScene()->FindChild("CaloOrthoDet1")->SetRnrState(kFALSE);

    // ... Turn ON rendering of detector in RPhi and RZ views
    calo2Dproj->fDetXYScene->FindChild("CaloOrthoDet0 [P]")->SetRnrState(kTRUE);
    calo2Dproj->fDetRZScene->FindChild("CaloOrthoDet1 [P]")->SetRnrState(kTRUE);

  }

  /*------------Function to create Tracker 2D tab:-------------*/
  void TEveMu2eMainWindow::PrepareTrackerProjectionTab(const art::Run& run){
    tracker2Dproj->fDetXYScene->DestroyElements();
    tracker2Dproj->fDetRZScene->DestroyElements();
    //fdetXY->DestroyElements();

    TEveElementList *orthodetXY = new TEveElementList("OrthoDetXY");
    TEveElementList *orthodetXZ = new TEveElementList("OrthoDetXZ");
    //TEveElementList *orthodetsplit = new TEveElementList("OrthoDet");

    TGeoVolume* topvol = geom->GetTopVolume();
    Mu2eTracker->DrawTrackerDetector(topvol, orthodetXZ, orthodetXY);
    //Mu2eTracker->DrawTrackerDetector(run, topvol, orthodetsplit);

    gEve->AddGlobalElement(orthodetXY);
    gEve->AddGlobalElement(orthodetXZ);

    // ... Import elements of the list into the projected views
    //TfXYMgr->ImportElements(orthodetsplit);
    //TfRZMgr->ImportElements(orthodetsplit);

    tracker2Dproj->fXYMgr->ImportElements(orthodetXY, tracker2Dproj->fDetXYScene);
    tracker2Dproj->fRZMgr->ImportElements(orthodetXZ, tracker2Dproj->fDetRZScene);

    // ... Turn OFF rendering of duplicate detector in main 3D view
    gEve->GetGlobalScene()->FindChild("OrthoDetXY")->SetRnrState(kFALSE);
    gEve->GetGlobalScene()->FindChild("OrthoDetXZ")->SetRnrState(kFALSE);
    //fdetXY->FindChild("OrthoDet [P]")->SetRnrState(kTRUE);
    // ... Turn ON rendering of detector in RPhi and RZ views
    tracker2Dproj->fDetXYScene->FindChild("OrthoDetXY [P]")->SetRnrState(kTRUE);
    tracker2Dproj->fDetRZScene->FindChild("OrthoDetXZ [P]")->SetRnrState(kTRUE);

  }

  /*------------Function to create CRV tab:-------------*/
  void TEveMu2eMainWindow::PrepareCRVProjectionTab(const art::Run& run){

  CRV2Dproj->fDetXYScene->DestroyElements();
  CRV2Dproj->fEvtXYScene->DestroyElements();

  TGeoVolume* topvol = geom->GetTopVolume();

  //TEveElementList *orthodetT1 = new TEveElementList("CRVT1OrthoDet");
  // TEveElementList *orthodetT2 = new TEveElementList("CRVT2OrthoDet");
  TEveElementList *orthodetT3 = new TEveElementList("CRVT3OrthoDet");
  TEveElementList *orthodetT4 = new TEveElementList("CRVT4OrthoDet");
  //TEveElementList *orthodetlist[] = {orthodetT1, orthodetT2, orthodetT3, orthodetT4};

  Mu2eCRV->DrawCRVDetector(run, topvol, orthodetT3, orthodetT4);

  //for (unsigned int i=0; i<2; i++){
  gEve->AddGlobalElement(orthodetT3);
  gEve->AddGlobalElement(orthodetT4);
  // }

  //for (unsigned int i=0; i<2; i++){
  CRV2Dproj->fXYMgr->ImportElements(orthodetT4, CRV2Dproj->fDetXYScene);
  CRV2Dproj->fRZMgr->ImportElements(orthodetT4, CRV2Dproj->fDetRZScene);
 // }

  // ... Turn OFF rendering of duplicate detector in main 3D view
  gEve->GetGlobalScene()->FindChild("CRVT4OrthoDet")->SetRnrState(kFALSE);
  gEve->GetGlobalScene()->FindChild("CRVT4OrthoDet")->SetRnrState(kFALSE);

  // ... Turn ON rendering of detector in RPhi and RZ views
  CRV2Dproj->fDetXYScene->FindChild("CRVT4OrthoDet [P]")->SetRnrState(kTRUE);
  CRV2Dproj->fDetRZScene->FindChild("CRVT4OrthoDet [P]")->SetRnrState(kTRUE);

  }

  /*------------Function to import the GDML and make 3D geometry:-------------*/
  void TEveMu2eMainWindow::SetRunGeometry(const art::Run& run, std::string gdmlname, int _diagLevel){
    if(gGeoManager){
      gGeoManager->GetListOfNodes()->Delete();
      gGeoManager->GetListOfVolumes()->Delete();
      gGeoManager->GetListOfShapes()->Delete();
    }
    gEve->GetGlobalScene()->DestroyElements();

    // Import the GDML of entire Mu2e Geometry
    ConfigFileLookupPolicy configFile;
    std::string fn =  configFile(gdmlname.c_str());

    geom = mu2e_geom->Geom_Interface::getGeom(fn);

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

    if(!GeomOpts.showbuilding){
      mu2e_geom->SolenoidsOnly(topnode);
      mu2e_geom->hideTop(topnode, _diagLevel);
    }
    if(GeomOpts.showDSOnly) mu2e_geom->InsideDS(topnode, false );
    if(GeomOpts.showInsidePS) mu2e_geom->InsidePS(topnode, false );
    if(GeomOpts.showCRV) mu2e_geom->InsideCRV(topnode, true);

    //Add static detector geometry to global scene
    gEve->AddGlobalElement(etopnode);
    geom->Draw("ogl");
  }

  /*------------Function to allow user to recheck the check box, and redraw data:-------------*/
  void TEveMu2eMainWindow::RedrawDataProducts(std::string type){
    if (type == "Clusters"){
      *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _emptydata.clustercol, calo2Dproj,true, fclustmin, fclustmax, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
    }
    if (type == "Hits"){
      if (_data.chcol !=0){*hitenergy = pass_data->AddComboHits(_firstLoop, _emptydata.chcol, tracker2Dproj,  true, fhitmin, fhitmax,ftimemin, ftimemax,_accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
      if (_data.tccol !=0){pass_data->AddTimeClusters(_firstLoop, _emptydata.tccol, tracker2Dproj,  true, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
      if (_data.chcol !=0){pass_data->AddTrkHits(_firstLoop, _emptydata.chcol, _emptydata.track_tuple, tracker2Dproj,  true, fhitmin, fhitmax,ftimemin, ftimemax,_accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
      if(_data.cryHitcol !=0){pass_data->AddCrystalHits(_firstLoop, _emptydata.cryHitcol, calo2Dproj, ftimemin, ftimemax, true, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);}

    }
    if (type == "Tracks"){
      pass_data->AddHelixPieceWise3D(_firstLoop, _emptydata.track_tuple, tracker2Dproj, ftimemin, ftimemax, true,  _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
    }
    if (type == "Cosmics"){
      if(_data.crvcoincol!= 0){pass_data->AddCRVInfo(_firstLoop, _emptydata.crvcoincol, ftimemin, ftimemax, CRV2Dproj, true,  _accumulate, TfXYMgr, TfRZMgr, proj4, proj5); }
      }
    if (type == "Cosmic Tracks"){std::cout<<"Cosmic tracks "<<std::endl;
      if(_data.cosmiccol!=0){pass_data->AddCosmicTrack(_firstLoop, _emptydata.cosmiccol, tracker2Dproj, ftimemin, ftimemax, true, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
    }
    if (type == "MC Trajectories"){
      if(_data.mctrajcol!=0){pass_mc->AddFullMCTrajectory(_firstLoop, _emptydata.mctrajcol, tracker2Dproj, true, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3, particles);}
    }
    gSystem->ProcessEvents();
    gClient->NeedRedraw(fTeRun);
    gApplication->Run(true);
    }

  /*------------Function to call button options:-------------*/
  Bool_t TEveMu2eMainWindow::ProcessMessage(Long_t msg, Long_t param1, Long_t param2){
    switch (GET_MSG(msg))
    {

    case kC_TEXTENTRY:
    switch (GET_SUBMSG(msg)){
      case kTE_TEXTCHANGED:
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
        //if (fhitmin < fhitmax) {*hitenergy = pass_data->AddComboHits(_firstLoop, _data.chcol, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
        //if (fhitmin < fhitmax) {pass_data->AddTrkHits(_firstLoop, _data.chcol, _data.track_tuple, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);}
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
            pass_data->AddTrkHits(_firstLoop, _data.chcol, _data.track_tuple, tracker2Dproj, false, fhitmin, fhitmax, ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
           }
           if(_data.tccol!=0) {
             pass_data->AddTimeClusters(_firstLoop, _data.tccol, tracker2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
           }if(_data.clustercol!=0){
            *clusterenergy = pass_data->AddCaloClusters(_firstLoop, _data.clustercol, calo2Dproj,  false, fclustmin, fclustmax,ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);
            }if(_data.crvcoincol!=0) {
               pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, ftimemin, ftimemax, CRV2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj4, proj5);
            }
          }
            }
      break;
    }
    break;
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
        pass_data->AddTimeClusters(_firstLoop, _data.tccol, tracker2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
        pass_data->AddTrkHits(_firstLoop, _data.chcol, _data.track_tuple, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
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
        pass_data->AddCRVInfo(_firstLoop, _data.crvcoincol, ftimemin, ftimemax, CRV2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj4, proj5);

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
          pass_mc->AddFullMCTrajectory(_firstLoop, _data.mctrajcol, tracker2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3, particles);
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

  /*------------Function to add the actual data to the plot (entrance from module):-------------*/
  void TEveMu2eMainWindow::setEvent(const art::Event& event, bool firstLoop, Data_Collections &data, double time, bool accumulate, int &runn, int &eventn, bool &eventSelected, bool isMCOnly)
  {

      _event=event.id().event();
      _subrun=event.id().subRun();
      _run=event.id().run();
      _firstLoop = firstLoop;
      _accumulate = accumulate;
      _data.chcol = data.chcol;
      _data.tccol = data.tccol;
      _data.clustercol = data.clustercol;
      _data.crvcoincol = data.crvcoincol;
      _data.track_tuple = data.track_tuple;
      _data.mctrajcol = data.mctrajcol;
      _data.crvcoincol = data.crvcoincol;
      _data.cryHitcol = data.cryHitcol;
      _data.cosmiccol = data.cosmiccol;

      std::string eveinfo = Form("Event : %i     Run : %i     Subrun : %i",_event,_run,_subrun);
      auto evinfo = new TEveText(eveinfo.c_str());
      double posy = -140.0;
      double posz = 0.0;
      evinfo->SetFontSize(15);
      evinfo->SetMainColor(kRed);
      evinfo->RefMainTrans().SetPos(posz,posy,posz);
      gEve->AddElement(evinfo);
    if(!isMCOnly){
      std::vector<const KalSeedCollection*> track_list = std::get<1>(data.track_tuple);
      std::vector<double> times = pass_data->getTimeRange(firstLoop, data.chcol, data.crvcoincol, data.clustercol, data.cryHitcol, DrawOpts.addCRVInfo, DrawOpts.addComboHits, DrawOpts.addClusters );

      if(DrawOpts.addCRVInfo){
        pass_data->AddCRVInfo(firstLoop, data.crvcoincol, ftimemin, ftimemax, CRV2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj4, proj5);
      }
      hitenergy = new vector<double>(2);

      if(DrawOpts.addComboHits) *hitenergy = pass_data->AddComboHits(firstLoop, data.chcol, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
      if(DrawOpts.addTrkHits) pass_data->AddTrkHits(firstLoop, data.chcol, data.track_tuple, tracker2Dproj, false, fhitmin, fhitmax,ftimemin, ftimemax, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
      if(DrawOpts.addTimeClusters) pass_data->AddTimeClusters(firstLoop, data.tccol, tracker2Dproj, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);

      clusterenergy = new std::vector<double>(2);

      if(DrawOpts.addClusters ) *clusterenergy = pass_data->AddCaloClusters(firstLoop, data.clustercol, calo2Dproj, false, fclustmin, fclustmax, ftimemin, ftimemax, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);

      if (DrawOpts.addCryHits) pass_data->AddCrystalHits(firstLoop, data.cryHitcol, calo2Dproj, ftimemin, ftimemax, false, _accumulate, CfXYMgr, CfRZMgr, proj0, proj1);

      if (DrawOpts.addTracks) pass_data->AddHelixPieceWise3D(firstLoop, data.track_tuple, tracker2Dproj,  ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);
      if (DrawOpts.addTracks) pass_data->FillKinKalTrajectory(firstLoop, data.track_tuple, tracker2Dproj,  TfXYMgr, TfRZMgr, proj2, proj3);
      

      if(DrawOpts.addCosmicTracks) pass_data->AddCosmicTrack(firstLoop, data.cosmiccol, tracker2Dproj, ftimemin, ftimemax, false, _accumulate, TfXYMgr, TfRZMgr, proj2, proj3);

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
    }

    if(DrawOpts.addMCTraj) pass_mc->AddFullMCTrajectory(firstLoop, data.mctrajcol, tracker2Dproj, false, _accumulate,  TfXYMgr, TfRZMgr, proj2, proj3, particles);

    gSystem->ProcessEvents();
    gSystem->IgnoreInterrupt();
    gSystem->IgnoreSignal(kSigTermination);
    gSystem->IgnoreSignal(kSigSegmentationViolation);

    gClient->NeedRedraw(fTeRun);

    gApplication->Run(true);

    gEve->Redraw3D(kTRUE);
    if(usereventSelected == true){
      eventn = eventToFind;
      runn = runToFind;
      eventSelected = true;
    }
    delete evinfo;
  }

  /*------------Function to find event:-------------*/
  int TEveMu2eMainWindow::getEventToFind(bool &findEvent) const
  {
    findEvent=_findEvent;
    return _eventToFind;
  }

  /*------------Function to check if closed:-------------*/
  bool TEveMu2eMainWindow::isClosed() const
  {
    return _isClosed;
  }
}

