//
// This class is inherited from ROOT's TPad class, and is used for _mainframe in EventDisplayFrame. It is used to catch updates to the display (e.g. rotations via mouse, zoom via TView3D context menu, etc.) so that the fields (phi, theta, xmin, xmax, ...) in the event display can get updated.
//
// $Id: EventDisplayPad.h,v 1.1 2011/02/05 01:20:08 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/05 01:20:08 $
//
// Original author Ralf Ehrlich
//

#ifndef EVENTDISPLAYPAD_H
#define EVENTDISPLAYPAD_H

#include <TPad.h>
#include <TView.h>
#include <TList.h>
#include "../EventDisplayFrame.h"

namespace mu2e_eventdisplay
{

class EventDisplayPad : public TPad
{
  EventDisplayPad(const EventDisplayPad &);
  EventDisplayPad& operator=(const EventDisplayPad &);

  EventDisplayFrame *_frame;

  public:
  EventDisplayPad() : TPad() 
  {
    _frame=NULL;
  }

  EventDisplayPad(const char* name, const char* title, Double_t xlow, Double_t ylow, Double_t xup, Double_t yup, Color_t color = -1, Short_t bordersize = -1, Short_t bordermode = -2) : TPad(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode) 
  {
    _frame=NULL;
  }

  void setEventDisplayFrame(EventDisplayFrame *frame) 
  {
    _frame=frame;
  }

  virtual void AbsCoordinates(Bool_t set)
  {
    if(_frame) _frame->fillZoomAngleFields();
    TPad::AbsCoordinates(set);
  }

  virtual ~EventDisplayPad() {}

  ClassDef(EventDisplayPad,0);
};

}
#endif
