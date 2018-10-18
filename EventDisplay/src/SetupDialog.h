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

  TGCheckButton *_checkButton[8];

  SetupDialog();
  SetupDialog(const SetupDialog &);
  SetupDialog& operator=(const SetupDialog &);

  public:
  SetupDialog(const TGWindow* p, EventDisplayFrame *eventDisplayFrame,
              bool whiteBackground, bool useHitColors, bool useTrackColors,
              bool showSupportStructures,
              bool showCRV,
              bool showOtherStructures,
              bool showMuonBeamStop, 
              bool showProtonAbsorber) : 
              TGMainFrame(p, 400, 400), _eventDisplayFrame(eventDisplayFrame)
  {
    SetCleanup(kDeepCleanup);

    _checkButton[0] = new TGCheckButton(this,"White background",11);
    _checkButton[1] = new TGCheckButton(this,"Use hit colors",14);
    _checkButton[2] = new TGCheckButton(this,"Use track colors",15);
    _checkButton[3] = new TGCheckButton(this,"Show support structures, calo disks, target",32);
    _checkButton[4] = new TGCheckButton(this,"Show CRV scintillator bars",36);
    _checkButton[5] = new TGCheckButton(this,"Show toy DS",37);
    _checkButton[6] = new TGCheckButton(this,"Show toy MBS",38);
    _checkButton[7] = new TGCheckButton(this,"Show MECO style (conical) proton absorber",39);

    TGLayoutHints *lh = new TGLayoutHints(kLHintsTop,2,1,2,2);
    for(int i=0; i<8; i++)
    {
      AddFrame(_checkButton[i], lh);
      _checkButton[i]->Associate(this);
    }

    _checkButton[0]->SetState(whiteBackground?kButtonDown:kButtonUp);
    _checkButton[1]->SetState(useHitColors?kButtonDown:kButtonUp);
    _checkButton[2]->SetState(useTrackColors?kButtonDown:kButtonUp);
    _checkButton[3]->SetState(showSupportStructures?kButtonDown:kButtonUp);
    _checkButton[4]->SetState(showCRV?kButtonDown:kButtonUp);
    _checkButton[5]->SetState(showOtherStructures?kButtonDown:kButtonUp);
    _checkButton[6]->SetState(showMuonBeamStop?kButtonDown:kButtonUp);
    _checkButton[7]->SetState(showProtonAbsorber?kButtonDown:kButtonUp);

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
                                            bool whiteBackground       =(_checkButton[0]->GetState()==kButtonDown);
                                            bool useHitColors          =(_checkButton[1]->GetState()==kButtonDown);
                                            bool useTrackColors        =(_checkButton[2]->GetState()==kButtonDown);
                                            bool showSupportStructures =(_checkButton[3]->GetState()==kButtonDown);
                                            bool showCRV               =(_checkButton[4]->GetState()==kButtonDown);
                                            bool showOtherStructures   =(_checkButton[5]->GetState()==kButtonDown);
                                            bool showMuonBeamStop      =(_checkButton[6]->GetState()==kButtonDown);
                                            bool showProtonAbsorber    =(_checkButton[7]->GetState()==kButtonDown);
                                            if(_eventDisplayFrame) _eventDisplayFrame->changeSetup(whiteBackground, 
                                                                                                   useHitColors, 
                                                                                                   useTrackColors,
                                                                                                   showSupportStructures,
                                                                                                   showCRV,
                                                                                                   showOtherStructures,
                                                                                                   showMuonBeamStop,
                                                                                                   showProtonAbsorber);
                                            CloseWindow();
                                          }
                                          break;
                       }
                       break;
    }
    return kTRUE;
  }

  ClassDef(SetupDialog,1);
};

}
#endif /* EventDisplay_src_SetupDialog_h */

