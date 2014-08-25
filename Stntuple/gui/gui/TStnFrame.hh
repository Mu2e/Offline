#ifndef TStnFrame_hh
#define TStnFrame_hh

#include "TGFrame.h"
#include "TCanvasImp.h"
#include "TRootEmbeddedCanvas.h"

class TControlBar;
class TGMenuBar;
class TGPopupMenu;
class TGLayoutHints;
class TGStatusBar;

class TStnVisManager;

//_____________________________________________________________________________
class TStnFrame: public TGMainFrame, public TCanvasImp {

protected:

  TGMenuBar           *fMenuBar;	// !
  TGPopupMenu         *fMenuFile;	// !
  TGPopupMenu         *fMenuEdit;	// !
  TGPopupMenu         *fMenuOption;	// !
  TGPopupMenu         *fMenuHelp;	// !

  TGLayoutHints       *fMenuBarLayout;	// !
  TGLayoutHints       *fMenuBarItemLayout; // !
  TGLayoutHints       *fMenuBarHelpLayout; // !


  TControlBar         *fEditorBar;       //! Editor control bar

  TRootEmbeddedCanvas *fEmbeddedCanvas;	// ! canvas widget
  TGLayoutHints       *fCanvasLayout;	// ! layout for canvas widget

  TGStatusBar*        fStatusBar;       // !
  TGLayoutHints       *fStatusBarLayout; // ! layout for the status bar

  TStnVisManager*     fVisManager;	// !
  int                 fView;		// !

public:
  TStnFrame(const char* name, const char* title, Int_t View, 
	    UInt_t w, UInt_t h,
	    UInt_t options = kMainFrame | kVerticalFrame);

  virtual ~TStnFrame();

					// ****** accessors

  TStnVisManager* GetVisManager() { return fVisManager; }
  int             GetView      () { return fView; }

					// ****** setters

  void SetVisManager(TStnVisManager* vm) { fVisManager = vm; }

					// ****** slots
  void DoOK();
  void DoCancel();
  void DoTab(Int_t id);
  void HandleButtons(Int_t id = -1);
  void EditorBar();
  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  void  ShowStatusBar(Bool_t show);
  void  SetStatusText(const char* txt = 0, Int_t partidx = 0);

  TCanvas*  GetCanvas() { return fEmbeddedCanvas->GetCanvas(); }

  virtual void CloseWindow();

  ClassDef(TStnFrame,1)
};


#endif
