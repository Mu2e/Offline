//
#ifndef Stntuple_gui_TStnVisFrame_hh
#ifndef Stntuple_gui_TStnVisFrame_hh

#include "TGMainFrame.h"

//_____________________________________________________________________________
class TMyMainFrame: public TGMainFrame {
public:
  TMyMainFrame();
  TMyMainFrame(const TGWindow* p, UInt_t w, UInt_t h, Int_t options);
  virtual ~TMyMainFrame() {};

  virtual Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
};

#endif
