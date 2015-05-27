//
#ifndef Stntuple_gui_TEvdMainFrame_hh
#define Stntuple_gui_TEvdMainFrame_hh

#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGStatusBar.h"
#include "Stntuple/gui/TStnVisManager.hh"


//_____________________________________________________________________________
class TEvdMainFrame: public TGMainFrame {
public:
  TEvdMainFrame();
  TEvdMainFrame(const TGWindow* p, UInt_t w, UInt_t h, Int_t options);
  virtual ~TEvdMainFrame() {};

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);

  TStnVisManager* vm;

  ClassDef(TEvdMainFrame, 0)
};

#endif
