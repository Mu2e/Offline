//-----------------------------------------------------------------------------
//  Mu2e Event display - stolen from TGeant visualization
//-----------------------------------------------------------------------------
#include <time.h>

#include "TROOT.h"
#include "TApplication.h"
#include "TVirtualX.h"

#include "TPad.h"
#include "TText.h"
#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGFileDialog.h"
#include "TControlBar.h"
#include "TInterpreter.h"
#include "TGStatusBar.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"

#include "Stntuple/gui/TStnFrame.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(TStnFrame)

//-----------------------------------------------------------------------------
enum TGeantCommandIdentifiers {
  M_FILE_OPEN,
  M_FILE_SAVE,
  M_FILE_SAVEAS,
  M_FILE_EXIT,

  M_TEST_DLG,
  M_TEST_MSGBOX,
  M_TEST_SLIDER,
  M_TEST_SHUTTER,
  M_TEST_PROGRESS,

  M_EDIT_EDITOR,
  M_EDIT_UNDO,
  M_EDIT_CLEARPAD,
  M_EDIT_CLEARCANVAS,

  M_OPTION_EVENT_STATUS,
  M_OPTION_AUTO_EXEC,
  M_OPTION_AUTO_RESIZE,
  M_OPTION_RESIZE_CANVAS,
  M_OPTION_MOVE_OPAQUE,
  M_OPTION_RESIZE_OPAQUE,
  M_OPTION_REFRESH,
  M_OPTION_STATISTICS,
  M_OPTION_HIST_TITLE,
  M_OPTION_FIT_PARAMS,
  M_OPTION_CAN_EDIT,

  M_HELP_CONTENTS,
  M_HELP_SEARCH,
  M_HELP_ABOUT
};

//-----------------------------------------------------------------------------
Int_t mb_button_id[9] = { kMBYes, kMBNo, kMBOk, kMBApply,
                          kMBRetry, kMBIgnore, kMBCancel,
                          kMBClose, kMBDismiss };

EMsgBoxIcon mb_icon[4] = { kMBIconStop, kMBIconQuestion,
                           kMBIconExclamation, kMBIconAsterisk };

const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            0,               0 };



//_____________________________________________________________________________
TStnFrame::TStnFrame(const char* Name,
		     const char* Title, 
		     Int_t       View,
		     UInt_t      w,
		     UInt_t      h, 
		     UInt_t      Options):

  TGMainFrame(gClient->GetRoot(),w, h, Options),
  fView(View)
{
  // create a TGeant event display window

//-----------------------------------------------------------------------------
//  create menu bar
//-----------------------------------------------------------------------------
  fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
				     0, 0, 1, 1);
  fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
  fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

  fMenuFile = new TGPopupMenu(gClient->GetRoot());
  fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
  fMenuFile->AddEntry("&Save", M_FILE_SAVE);
  fMenuFile->AddEntry("S&ave as...", M_FILE_SAVEAS);
  fMenuFile->AddEntry("&Close", -1);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("&Print", -1);
  fMenuFile->AddEntry("P&rint setup...", -1);
  fMenuFile->AddSeparator();
  fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

  fMenuFile->DisableEntry(M_FILE_SAVEAS);

  fMenuHelp = new TGPopupMenu(gClient->GetRoot());
  fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
  fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
  fMenuHelp->AddSeparator();
  fMenuHelp->AddEntry("&About", M_HELP_ABOUT);

					// Menu item EDIT

  fMenuEdit = new TGPopupMenu(gClient->GetRoot());
  fMenuEdit->AddEntry("&Editor",             M_EDIT_EDITOR);
  fMenuEdit->AddEntry("&Undo",               M_EDIT_UNDO);
  fMenuEdit->AddEntry("Clear &Pad",          M_EDIT_CLEARPAD);
  fMenuEdit->AddEntry("&Clear Canvas",       M_EDIT_CLEARCANVAS);

  fMenuOption = new TGPopupMenu(gClient->GetRoot());
  fMenuOption->AddEntry("&Event Status",         M_OPTION_EVENT_STATUS);
  fMenuOption->AddEntry("&Pad Auto Exec",        M_OPTION_AUTO_EXEC);
  fMenuOption->AddSeparator();
  fMenuOption->AddEntry("&Auto Resize Canvas",   M_OPTION_AUTO_RESIZE);
  fMenuOption->AddEntry("&Resize Canvas",        M_OPTION_RESIZE_CANVAS);
  fMenuOption->AddEntry("&Move Opaque",          M_OPTION_MOVE_OPAQUE);
  fMenuOption->AddEntry("Resize &Opaque",        M_OPTION_RESIZE_OPAQUE);
  fMenuOption->AddEntry("R&efresh",              M_OPTION_REFRESH);
  fMenuOption->AddSeparator();
  fMenuOption->AddEntry("Show &Statistics",      M_OPTION_STATISTICS);
  fMenuOption->AddEntry("Show &Histogram Title", M_OPTION_HIST_TITLE);
  fMenuOption->AddEntry("Show &Fit Parameters",  M_OPTION_FIT_PARAMS);
  fMenuOption->AddEntry("Can Edit Histograms",   M_OPTION_CAN_EDIT);

					// define menu handlers
  fMenuFile->Associate(this);
  fMenuEdit->Associate(this);
  fMenuOption->Associate(this);
  fMenuHelp->Associate(this);

  fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);

  fMenuBar->AddPopup("&File"  , fMenuFile  , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Edit"  , fMenuEdit  , fMenuBarItemLayout);
  fMenuBar->AddPopup("&Option", fMenuOption, fMenuBarItemLayout);
  fMenuBar->AddPopup("&Help"  , fMenuHelp  , fMenuBarHelpLayout);
  
  AddFrame(fMenuBar, fMenuBarLayout);
//-----------------------------------------------------------------------------
// Create canvas and canvas container that will host the ROOT graphics
//-----------------------------------------------------------------------------
  fEmbeddedCanvas = new TRootEmbeddedCanvas(Name, 
					    this, 
					    GetWidth()+4, 
					    GetHeight()+4,
					    kSunkenFrame | kDoubleBorder);

  fCanvasLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
  TCanvas* c = fEmbeddedCanvas->GetCanvas();

					// pad for the event display
  c->cd();
  TString name1(Name);
  name1 += "_1";
  TPad *p1 = new TPad(name1, "p1",0.0,0.0,1,0.96);
  p1->Draw();
  c->cd();

  TString name2(Name);
  name2 += "_2";
  TPad *p2 = new TPad(name2, "p2",0.7,0.96,1,1);
  p2->Draw();
  p2->cd();

  time_t t2 = time(0);
  tm* t22 = localtime(&t2);
  TText *text = new TText(0.05,0.3,asctime(t22));
  text->SetTextFont(22);
  text->SetTextSize(0.4);
  text->Draw();
  p2->Modified();

  c->cd();

  TString name_title(Name);
  name2 += "_title";
  TPad *title = new TPad(name_title, "title",0.,0.96,0.7,1);
  title->Draw();

  title->Modified();
 

  AddFrame(fEmbeddedCanvas, fCanvasLayout);
//-----------------------------------------------------------------------------
// no editor bar by default
//-----------------------------------------------------------------------------
  fEditorBar     = NULL;
//-----------------------------------------------------------------------------
// status bar (hidden by default)
//-----------------------------------------------------------------------------
   int parts[] = { 33, 10, 10, 47 };
   fStatusBar = new TGStatusBar(this, 10, 10);
   fStatusBar->SetParts(parts, 4);
   fStatusBarLayout = new TGLayoutHints(kLHintsBottom | kLHintsLeft | 
					kLHintsExpandX, 2, 2, 1, 1);
   AddFrame(fStatusBar, fStatusBarLayout);
   HideFrame(fStatusBar);
//-----------------------------------------------------------------------------
// final actions
//-----------------------------------------------------------------------------
  MapSubwindows();
  Resize(GetDefaultSize());

  SetWindowName(Title);

  MapWindow();
}


//_____________________________________________________________________________
TStnFrame::~TStnFrame() {
  // Delete window, as TStnFrame's do not exist by themselves, but they are
  // always managed by TVisManager, we need to erase this frame from the
  // list of frames

  TVisManager::Instance()->GetListOfCanvases()->Remove(GetCanvas());

  delete fStatusBar;
  delete fStatusBarLayout;

  delete fEmbeddedCanvas;
  delete fCanvasLayout;

  delete fMenuFile;
  delete fMenuEdit;
  delete fMenuOption;
  delete fMenuHelp;

  delete fMenuBarItemLayout;
  delete fMenuBarHelpLayout;
  delete fMenuBarLayout;

  delete fMenuBar;

  if (fEditorBar) delete fEditorBar;

}

//______________________________________________________________________________
void TStnFrame::EditorBar() {
  // Create the Editor Controlbar

   TControlBar *ed = new TControlBar("vertical", "Editor");
   ed->AddButton("Arc",       "gROOT->SetEditorMode(\"Arc\")",       "Create an arc of circle");
   ed->AddButton("Line",      "gROOT->SetEditorMode(\"Line\")",      "Create a line segment");
   ed->AddButton("Arrow",     "gROOT->SetEditorMode(\"Arrow\")",     "Create an Arrow");
   ed->AddButton("Button",    "gROOT->SetEditorMode(\"Button\")",    "Create a user interface Button");
   ed->AddButton("Diamond",   "gROOT->SetEditorMode(\"Diamond\")",   "Create a diamond");
   ed->AddButton("Ellipse",   "gROOT->SetEditorMode(\"Ellipse\")",   "Create an Ellipse");
   ed->AddButton("Pad",       "gROOT->SetEditorMode(\"Pad\")",       "Create a pad");
   ed->AddButton("Pave",      "gROOT->SetEditorMode(\"Pave\")",      "Create a Pave");
   ed->AddButton("PaveLabel", "gROOT->SetEditorMode(\"PaveLabel\")", "Create a PaveLabel (prompt for label)");
   ed->AddButton("PaveText",  "gROOT->SetEditorMode(\"PaveText\")",  "Create a PaveText");
   ed->AddButton("PavesText", "gROOT->SetEditorMode(\"PavesText\")", "Create a PavesText");
   ed->AddButton("PolyLine",  "gROOT->SetEditorMode(\"PolyLine\")",  "Create a PolyLine (TGraph)");
   ed->AddButton("CurlyLine", "gROOT->SetEditorMode(\"CurlyLine\")", "Create a Curly/WavyLine");
   ed->AddButton("CurlyArc",  "gROOT->SetEditorMode(\"CurlyArc\")",  "Create a Curly/WavyArc");
   ed->AddButton("Text/Latex","gROOT->SetEditorMode(\"Text\")",      "Create a Text/Latex string");
   ed->AddButton("Marker",    "gROOT->SetEditorMode(\"Marker\")",    "Create a marker");
   ed->Show();
   fEditorBar = ed;
}

//_____________________________________________________________________________
Bool_t TStnFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) {
   // Handle menu items.

  TCanvas* c;
  int message = GET_MSG(msg);
//   double     x,y;
//   int        px, py;

  TVisManager* vm = TVisManager::Instance();

  c = GetCanvas();

  switch (message) {
  case kC_COMMAND:
    int submessage = GET_SUBMSG(msg);

    printf(" *** TStnFrame::ProcessMessage SUBMESSAGE: %i \n",submessage);

    switch (submessage) {
    case kCM_MENU:
      switch (parm1) {
 //-----------------------------------------------------------------------------
//  FILE menu
//-----------------------------------------------------------------------------
      case M_FILE_OPEN: {
	TGFileInfo fi;
	fi.fFileTypes = (const char **)filetypes;
	new TGFileDialog(gClient->GetRoot(), this, kFDOpen,&fi);
      }
      break;

      case M_FILE_SAVE:
	printf("M_FILE_SAVE\n");
	break;
	
      case M_FILE_EXIT:
	CloseWindow();   // terminate theApp no need to use SendCloseMessage()
	break;
//-----------------------------------------------------------------------------
//  EDIT menu
//-----------------------------------------------------------------------------
      case M_EDIT_EDITOR:
	EditorBar();
	break;
//-----------------------------------------------------------------------------
//  OPTIONS menu
//-----------------------------------------------------------------------------
      case M_OPTION_EVENT_STATUS:

	printf(" *** TStnFrame::ProcessMessage M_OPTION_EVENT_STATUS: msg = %li parm1 = %li parm2 = %li\n", 
	       msg,parm1,parm2);
	c->ToggleEventStatus();
	if (c->GetShowEventStatus()) {
	  ShowFrame(fStatusBar);
	  fMenuOption->CheckEntry(M_OPTION_EVENT_STATUS);
	} 
	else {
	  HideFrame(fStatusBar);
	  fMenuOption->UnCheckEntry(M_OPTION_EVENT_STATUS);
	}
	break;

      case M_OPTION_CAN_EDIT:
	if (gROOT->GetEditHistograms()) {
	  gROOT->SetEditHistograms(kFALSE);
	  fMenuOption->UnCheckEntry(M_OPTION_CAN_EDIT);
	} 
	else {
	  gROOT->SetEditHistograms(kTRUE);
	  fMenuOption->CheckEntry(M_OPTION_CAN_EDIT);
	}
	break;
//-----------------------------------------------------------------------------
//  default
//-----------------------------------------------------------------------------
      default:
	if (vm->DebugLevel() > 0) {
	  printf(" *** TStnFrame::ProcessMessage msg = %li parm1 = %li parm2 = %li\n", 
		 msg,parm1,parm2);
	}
	break;
      }
    default:
      if (vm->DebugLevel() > 0) {
	printf(" *** TStnFrame::ProcessMessage msg = %li parm1 = %li parm2 = %li\n", 
	       msg,parm1,parm2);
      }
      break;
    }
  }

//   px = gPad->GetEventX();
//   py = gPad->GetEventY();
  
//   x = ((TPad*) gPad)->AbsPixeltoX(px);
//   y = ((TPad*) gPad)->AbsPixeltoY(py);
  
//   fStatusBar->SetText(Form("Z = %8.3f",x),0);
//   fStatusBar->SetText(Form("R = %8.3f",y),1);
//   fStatusBar->SetText("ccc",2);
//   fStatusBar->SetText("ddd",3);
  
  return true;
}



//_____________________________________________________________________________
void TStnFrame::CloseWindow() {
   // Called when window is closed via the window manager.

   TVirtualPad *savepad = gPad;
   gPad = 0;        // hide gPad from CINT
   gInterpreter->DeleteGlobal(fEmbeddedCanvas->GetCanvas());
   gPad = savepad;  // restore gPad for ROOT

   delete this;
}


//______________________________________________________________________________
void TStnFrame::ShowStatusBar(Bool_t show)
{
   // Show or hide statusbar.

   if (show) {
      ShowFrame(fStatusBar);
      fMenuOption->CheckEntry(M_OPTION_EVENT_STATUS);
   } else {
      HideFrame(fStatusBar);
      fMenuOption->UnCheckEntry(M_OPTION_EVENT_STATUS);
   }
}

//______________________________________________________________________________
void TStnFrame::SetStatusText(const char *txt, Int_t partidx)
{
   // Set text in status bar.

   fStatusBar->SetText(txt, partidx);
}

//_____________________________________________________________________________
void TStnFrame::DoOK()
{
   printf("\nTerminating dialog: OK pressed\n");

   // Send a close message to the main frame. This will trigger the
   // emission of a CloseWindow() signal, which will then call
   // TStnFrame::CloseWindow(). Calling directly CloseWindow() will cause
   // a segv since the OK button is still accessed after the DoOK() method.
   // This works since the close message is handled synchronous (via
   // message going to/from X server).

   SendCloseMessage();

   // The same effect can be obtained by using a singleshot timer:
   //TTimer::SingleShot(50, "TStnFrame", this, "CloseWindow()");
}

//_____________________________________________________________________________
void TStnFrame::DoCancel() {
  printf("\nTerminating dialog: Cancel pressed\n");
  SendCloseMessage();
}

//_____________________________________________________________________________
void TStnFrame::HandleButtons(Int_t id) {
  // Handle different buttons.
}

//_____________________________________________________________________________
void TStnFrame::DoTab(Int_t id) {
   printf("*** TStnFrame::DoTab : Tab item %d activated\n", id);
}

