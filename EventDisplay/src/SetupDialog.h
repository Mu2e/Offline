#ifndef EventDisplay_src_SetupDialog_h
#define EventDisplay_src_SetupDialog_h

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include "EventDisplay/src/EventDisplayFrame.h"

namespace mu2e_eventdisplay
{

class SetupDialog : public TGMainFrame
{
  EventDisplayFrame *_eventDisplayFrame;

  TGCheckButton *_checkButton1;
  TGCheckButton *_checkButton2;
  TGCheckButton *_checkButton3;
  TGCheckButton *_checkButton4;
  TGCheckButton *_checkButton5;

  SetupDialog();
  SetupDialog(const SetupDialog &);
  SetupDialog& operator=(const SetupDialog &);

  public:
  SetupDialog(const TGWindow* p, EventDisplayFrame *eventDisplayFrame,
               bool whiteBackground, bool showUnhitStraws, bool showUnhitCrystals,
               bool useHitColors, bool useTrackColors) :
               TGMainFrame(p, 400, 400), _eventDisplayFrame(eventDisplayFrame)
  {
    SetCleanup(kDeepCleanup);

    _checkButton1 = new TGCheckButton(this,"White background",11);
    _checkButton2 = new TGCheckButton(this,"Show unhit straws",12);
    _checkButton3 = new TGCheckButton(this,"Show unhit calorimeter crystals",13);
    _checkButton4 = new TGCheckButton(this,"Use hit colors",14);
    _checkButton5 = new TGCheckButton(this,"Use track colors",15);

    TGLayoutHints *lh = new TGLayoutHints(kLHintsTop,2,1,2,2);
    AddFrame(_checkButton1, lh);
    AddFrame(_checkButton2, lh);
    AddFrame(_checkButton3, lh);
    AddFrame(_checkButton4, lh);
    AddFrame(_checkButton5, lh);

    _checkButton1->Associate(this);
    _checkButton2->Associate(this);
    _checkButton3->Associate(this);
    _checkButton4->Associate(this);
    _checkButton5->Associate(this);

    _checkButton1->SetState(whiteBackground?kButtonDown:kButtonUp);
    _checkButton2->SetState(showUnhitStraws?kButtonDown:kButtonUp);
    _checkButton3->SetState(showUnhitCrystals?kButtonDown:kButtonUp);
    _checkButton4->SetState(useHitColors?kButtonDown:kButtonUp);
    _checkButton5->SetState(useTrackColors?kButtonDown:kButtonUp);

    TGHorizontalFrame *subFrame1  = new TGHorizontalFrame(this,500,20);
    AddFrame(subFrame1,  lh);
    TGTextButton *button1 = new TGTextButton(subFrame1, "&Ok", 101);
    TGTextButton *button2 = new TGTextButton(subFrame1, "&Cancel", 102);
    subFrame1->AddFrame(button1, lh);
    subFrame1->AddFrame(button2, lh);
    button1->Associate(this);
    button2->Associate(this);

    SetWindowName("Setup");
    MapSubwindows();
    MapWindow();
    Layout();

    gClient->WaitFor(this);
  }

  virtual ~SetupDialog()
  {
  }

  Bool_t ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
  {
    switch(GET_MSG(msg))
    {
      case kC_COMMAND: switch(GET_SUBMSG(msg))
                       {
                         case kCM_BUTTON: if(param1==102) CloseWindow();
                                          if(param1==101)
                                          {
                                            bool whiteBackground  =(_checkButton1->GetState()==kButtonDown);
                                            bool showUnhitStraws  =(_checkButton2->GetState()==kButtonDown);
                                            bool showUnhitCrystals=(_checkButton3->GetState()==kButtonDown);
                                            bool useHitColors     =(_checkButton4->GetState()==kButtonDown);
                                            bool useTrackColors   =(_checkButton5->GetState()==kButtonDown);
                                            if(_eventDisplayFrame) _eventDisplayFrame->changeSetup(whiteBackground, 
                                                                                                   showUnhitStraws, showUnhitCrystals,
                                                                                                   useHitColors, useTrackColors);
                                            CloseWindow();
                                          }
                                          break;
                       }
                       break;
    }
    return kTRUE;
  }
};

}
#endif /* EventDisplay_src_SetupDialog_h */

