#ifndef TStnVisManager_hh
#define TStnVisManager_hh

#include "TObjArray.h"
#include "Stntuple/base/TVisManager.hh"

#ifndef __CINT__

#include "art/Framework/Principal/Event.h"

#else

namespace art {
  class Event;
}

#endif

class TControlBar;
class TGMenuBar;
class TGPopupMenu;
class TGLayoutHints;
class TGMainFrame;

class TTrkXYView;
class TTrkRZView;
class TCalView;
class TSubdetector;
class TExtrapolator;

class TStnVisManager: public TVisManager {
public:
  enum { 
    kXYView       = 1,
    kRZView       = 2,
    kCalView      = 3
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TGMainFrame*        fMain;
  TGMenuBar           *fMenuBar;	  // !
  TGPopupMenu         *fMenuSubdetectors; // !
  TGPopupMenu         *fMenuHelp;	  // !

  TGLayoutHints       *fMenuBarLayout;	  // !
  TGLayoutHints       *fMenuBarItemLayout; // !
  TGLayoutHints       *fMenuBarHelpLayout; // !

					// vis. manager also holds a list of
					// objects to be displayed: has to be
					// the same for all the views

  TObjArray*          fListOfDetectors;
  TSubdetector*       fClosestSubdetector;

  TTrkXYView*         fTrkXYView;
  TTrkRZView*         fTrkRZView;
  TCalView*           fCalView[4];	// to provide for the now obsolete 
					// vane-based geometry

  TExtrapolator*      fExtrapolator;

  const art::Event*   fEvent;

  int                 fMinStation;
  int                 fMaxStation;
  int                 fTimePeak;
  int                 fDebugLevel;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TStnVisManager(const char* name  = "TStnVisManager",
		 const char* title = "TStnVisManager");

  virtual ~TStnVisManager();

  static TStnVisManager* Instance();
					// ****** accessors

  TSubdetector*  GetClosestSubdetector() { return fClosestSubdetector; }
  TExtrapolator* GetExtrapolator      () { return fExtrapolator; }

  TObjArray*     GetListOfDetectors() { return fListOfDetectors; }

  void          AddDetector(TObject* det) { fListOfDetectors->Add(det); }

  const art::Event* Event() { return fEvent; }

  int  MinStation() { return fMinStation; }
  int  MaxStation() { return fMaxStation; }
  int  TimePeak  () { return fTimePeak;   }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetEvent(art::Event& Evt) { fEvent = &Evt; }

  void SetClosestSubdetector(TSubdetector* det) { fClosestSubdetector = det; }
  void SetExtrapolator      (TExtrapolator*  x) { fExtrapolator       = x; }

  void SetStations(int IMin, int IMax);
  void SetTimePeak(int I);

  virtual TCanvas*  NewCanvas(const char* Name, 
			      const char* Title, 
			      Int_t       SizeX, 
			      Int_t       SizeY);

  Int_t   OpenTrkXYView();
  Int_t   OpenTrkXYView(TTrkXYView* mother,
			Axis_t x1, Axis_t y1, 
			Axis_t x2, Axis_t y2);
  
  Int_t   OpenTrkRZView();
  Int_t   OpenTrkRZView(TTrkRZView* mother,
			Axis_t x1, Axis_t y1, 
			Axis_t x2, Axis_t y2);
  
  Int_t   OpenCalView();
  Int_t   OpenCalView(TObject* Mother,
		      Axis_t x1, Axis_t y1, 
		      Axis_t x2, Axis_t y2);


  void    CloseWindow();

  ClassDef(TStnVisManager,0)
};
#endif
