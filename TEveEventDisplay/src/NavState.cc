//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>

#include "TROOT.h"
#include "TEveEventDisplay/src/dict_classes/NavState.h"
using namespace mu2e;

  static int gsNavState    = 0;
  static int gsTargetRun   = 0;
  static int gsTargetEvent = 0;

  //......................................................................

  int NavState::Which() { return gsNavState; }

  //......................................................................

  void NavState::Set(int which)
  {
    gsNavState = which;
    gROOT->GetApplication()->Terminate();
  }

  //......................................................................

  void NavState::SetTarget(int run, int event)
  {
    gsTargetRun = run;
    gsTargetEvent = event;
  }

  //......................................................................

  int NavState::TargetRun() { return gsTargetRun; }

  //......................................................................

  int NavState::TargetEvent() { return gsTargetEvent; }

