///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TVirtualX.h"

#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGStatusBar.h"

#include "Stntuple/gui/TTrkXYView.hh"
#include "Stntuple/gui/TTrkRZView.hh"
#include "Stntuple/gui/TCalView.hh"
#include "Stntuple/gui/TStnFrame.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(TStnVisManager)

//-----------------------------------------------------------------------------
enum TStmVisManagerCommandIdentifiers {
  M_TRACKER_XY,
  M_TRACKER_RZ,
  M_CALORIMETER_XY,
  M_EXIT,

  M_OPTION_EVENT_STATUS,

  M_HELP_CONTENTS,
  M_HELP_SEARCH,
  M_HELP_ABOUT
};


//_____________________________________________________________________________
class TMyMainFrame: public TGMainFrame {
public:
  TMyMainFrame();
  TMyMainFrame(const TGWindow* p, UInt_t w, UInt_t h, Int_t options);
  virtual ~TMyMainFrame() {};

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
};

//_____________________________________________________________________________
TMyMainFrame::TMyMainFrame():
  TGMainFrame(gClient->GetRoot(),100,200,kMainFrame)
{}

//_____________________________________________________________________________
TMyMainFrame::TMyMainFrame(const TGWindow* p, UInt_t w, UInt_t h, Int_t options):
  TGMainFrame(gClient->GetRoot(),w,h,options)
{}




//_____________________________________________________________________________
TStnVisManager::TStnVisManager(const char* Name, const char* Title):
  TVisManager(Name,Title)
{
  if (gROOT->IsBatch()) return ;

  fMain = new  TMyMainFrame(gClient->GetRoot(),200,100, 
			    kMainFrame | kVerticalFrame);
//-----------------------------------------------------------------------------
//  create menu bar
//-----------------------------------------------------------------------------
  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
				     0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fMenuSubdetectors = new TGPopupMenu(gClient->GetRoot());
  fMenuSubdetectors->AddEntry("&Tracker(XY view)", M_TRACKER_XY);
  fMenuSubdetectors->AddEntry("&Tracker(RZ view)", M_TRACKER_RZ);
  fMenuSubdetectors->AddEntry("&Calorimeter"     , M_CALORIMETER_XY);
  fMenuSubdetectors->AddEntry("&Exit"            , M_EXIT);
  fMenuSubdetectors->Associate(fMain);

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&Contents" , M_HELP_CONTENTS);
  fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("&About"    , M_HELP_ABOUT);
  fMenuHelp->Associate(fMain);

  fMenuBar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
  fMenuBar->AddPopup("&Subdetectors", fMenuSubdetectors, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help"        , fMenuHelp        , fMenuBarHelpLayout);

  fMain->AddFrame(fMenuBar, fMenuBarLayout);
//-----------------------------------------------------------------------------
// views
//-----------------------------------------------------------------------------
  fTrkXYView       = new TTrkXYView();
					// for vanes - 4 "views"
  fCalView[0]      = new TCalView (0);
  fCalView[1]      = new TCalView (1);

  fListOfDetectors = new TObjArray(10);
//-----------------------------------------------------------------------------
// final actions
//-----------------------------------------------------------------------------
  fMain->MapSubwindows();
  fMain->Resize(fMain->GetDefaultSize());
  fMain->Resize(200,100);

  fMain->SetWindowName(Title);

  fMain->MapWindow();

  fMinStation = 0;
  fMaxStation = 50;
  
  fTimePeak   = -1;
}



//_____________________________________________________________________________
TStnVisManager::~TStnVisManager() {

  if (! gROOT->IsBatch()) {

    delete fTrkXYView;
					// only two views for disk calorimeter
    delete fCalView[0];
    delete fCalView[1];
    
    delete fMenuBarHelpLayout;
    
    delete fMenuBarItemLayout;
    delete fMenuSubdetectors;
    
    delete fMenuBarLayout;
    delete fMenuBar;
    
    delete fMain;
    
    delete fListOfDetectors;
  }
}

//_____________________________________________________________________________
TStnVisManager* TStnVisManager::Instance() {
  if (TVisManager::fgInstance != NULL) {
    return (TStnVisManager*) TVisManager::fgInstance;
  }
  else {
    return new TStnVisManager();
  }
}

//_____________________________________________________________________________
TCanvas* TStnVisManager::NewCanvas(const char* Name, 
				   const char* Title, 
				   Int_t       SizeX, 
				   Int_t       SizeY)
{
  TStnFrame* win = new TStnFrame(Name,Title,0,SizeX,SizeY);
  TCanvas*c = win->GetCanvas();
  DeclareCanvas(c);
  return c;
}


//_____________________________________________________________________________
Bool_t TMyMainFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) {
   // Handle menu items.

  int message = GET_MSG(msg);

  switch (message) {
  case kC_COMMAND:
    int submessage = GET_SUBMSG(msg);
    switch (submessage) {
    case kCM_MENU:
      switch (parm1) {
//-----------------------------------------------------------------------------
//  SUBDETECTOR menu
//-----------------------------------------------------------------------------
      case M_TRACKER_XY: 
	TStnVisManager::Instance()->OpenTrkXYView();
	break;

      case M_CALORIMETER_XY: 
	TStnVisManager::Instance()->OpenCalView();
	break;

      case M_EXIT:
	TStnVisManager::Instance()->CloseWindow();
	break;
//-----------------------------------------------------------------------------
//  default
//-----------------------------------------------------------------------------
      default:
//  	printf(" *** TStnFrame::ProcessMessage msg = %i parm1 = %i parm2 = %i\n", 
//  	       msg,parm1,parm2);
	break;
      }
    default:
//      printf(" *** TStnFrame::ProcessMessage msg = %i parm1 = %i parm2 = %i\n", 
//     msg,parm1,parm2);
      break;
    }
  }
  return true;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenTrkXYView() {
   // open new XY view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name,"xy_view_%i",n);
  sprintf(title,"XY view number %i",n);
  
  TStnFrame* win = new TStnFrame(name, title, kXYView,740,760);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(-800.,-800.,800.,800.);
  p1->cd();
  fTrkXYView->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenTrkXYView(TTrkXYView* mother, Axis_t x1, Axis_t y1,
				    Axis_t x2, Axis_t y2) 
{
   // open new XY view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name,"xy_view_%i",n);
  sprintf(title,"XY view number %i",n);

				// try to preserve the aspect ratio
  Int_t   xsize, ysize;

  xsize = 540;
  ysize = (Int_t) (xsize*TMath::Abs((y2-y1)/(x2-x1))+20);

  TStnFrame* win = new TStnFrame(name, title, kXYView, xsize,ysize);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(x1,y1,x2,y2);
  p1->cd();
  mother->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//-----------------------------------------------------------------------------
// open new RZ view of the detector with the default options
//-----------------------------------------------------------------------------
Int_t TStnVisManager::OpenTrkRZView() {

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name ,"rz_view_%i"       ,n);
  sprintf(title,"RZ view number %i",n);
  
  TStnFrame* win = new TStnFrame(name, title, kRZView,740,760);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(-800.,-800.,800.,800.);
  p1->cd();
  fTrkRZView->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//-----------------------------------------------------------------------------
// open new RZ view of the detector with the default options
//-----------------------------------------------------------------------------
Int_t TStnVisManager::OpenTrkRZView(TTrkRZView* mother, 
				    Axis_t      x1, Axis_t y1,
				    Axis_t      x2, Axis_t y2) 
{

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name ,"rz_view_%i"       ,n);
  sprintf(title,"RZ view number %i",n);
//-----------------------------------------------------------------------------
// try to preserve the aspect ratio
//-----------------------------------------------------------------------------
  Int_t   xsize, ysize;

  xsize = 540;
  ysize = (Int_t) (xsize*TMath::Abs((y2-y1)/(x2-x1))+20);

  TStnFrame* win = new TStnFrame(name, title, kRZView, xsize,ysize);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(x1,y1,x2,y2);
  p1->cd();
  mother->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}


//_____________________________________________________________________________
Int_t TStnVisManager::OpenCalView() {
  // open new calorimeter view of the detector with the default options
  // start from the disk calorimeter

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name ,"cal_view_%i",n);
  sprintf(title,"CAL view number %i",n);
  
  TStnFrame* win = new TStnFrame(name,title, kCalView, 1200, 600);
  TCanvas*   c   = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
//-----------------------------------------------------------------------------
// the disk calorimeter view has two pads, one per disk
// the vane-based calorimeter display should have 4 pads in this view
//-----------------------------------------------------------------------------
					// divide horisontally
  p1->Divide(2,1);
					// ranges in mm
  p1->cd(1);
  gPad->Range(-800.,-800.,800.,800.);
  fCalView[0]->SetPad(gPad);
  fCalView[0]->Draw("cal");
  gPad->Modified();

  p1->cd(2);
  gPad->Range(-800.,-800.,800.,800.);
  fCalView[1]->SetPad(gPad);
  fCalView[1]->Draw("cal");
  gPad->Modified();
					// draw title
  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
   c->Update();
  return 0;
}

//_____________________________________________________________________________
Int_t TStnVisManager::OpenCalView(TObject* Mother, Axis_t x1, Axis_t y1,
				  Axis_t x2, Axis_t y2) 
{
   // open new calorimeter view of the detector with the default options

//   int n = fListOfCanvases->GetSize();

//   char name[100], title[100];

//   sprintf(name,"ces_view_%i",n);
//   sprintf(title,"CES view number %i",n);

// 				// try to preserve the aspect ratio
//   Int_t   xsize, ysize;

//   xsize = 540;
//   ysize = (Int_t) (xsize*TMath::Abs((y2-y1)/(x2-x1))+20);

//   TStnFrame* win = new TStnFrame(name, title, kCesStripView, xsize,ysize);
//   TCanvas* c = win->GetCanvas();
//   fListOfCanvases->Add(c);
//   c->Divide(1,2);

//   TString name1(name);
//   name1 += "_1";
//   TPad* p1 = (TPad*) c->FindObject(name1);

//   p1->Divide(2,1);

//   p1->cd(1);
//   gPad->Range(x1,y1,x2,y2);

//   fCalSectionView[0]->Draw();

//   p1->cd(2);
//   gPad->Range(x1,y1,x2,y2);
//   fCalSectionView[1]->Draw();

//   TString name_title(name);
//   name1 += "_title";
//   TPad* title_pad = (TPad*) c->FindObject(name_title);
//   title_pad->cd();
//   fTitleNode->Draw();

//   c->Modified();
//   c->Update();
  return 0;
}

//_____________________________________________________________________________
void TStnVisManager::CloseWindow() {
   // Called when window is closed via the window manager.

   delete this;
}


//-----------------------------------------------------------------------------
void TStnVisManager::SetStations(int IMin, int IMax) {
  fMinStation = IMin;
  fMaxStation = IMax;
}

//-----------------------------------------------------------------------------
void TStnVisManager::SetTimePeak(int I) {
  fTimePeak = I;
}

