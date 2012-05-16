#ifndef EventDisplay_src_FilterDialog_h
#define EventDisplay_src_FilterDialog_h

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include "DataInterface.h"

namespace mu2e_eventdisplay
{

class FilterDialog : public TGMainFrame
{
  DataInterface *_dataInterface;

  TGTextEntry   *_textEntry1;
  TGTextEntry   *_textEntry2;
  TGTextEntry   *_textEntry3;
  TGTextEntry   *_textEntry4;
  TGCheckButton *_checkButton1;
  TGCheckButton *_checkButton2;
  TGCheckButton *_checkButton3;
  TGCheckButton *_checkButton4;
  TGCheckButton *_checkButton5;
  TGCheckButton *_checkButton6;

  FilterDialog();
  FilterDialog(const FilterDialog &);
  FilterDialog& operator=(const FilterDialog &);

  public:
  FilterDialog(const TGWindow* p, DataInterface *dataInterface,
               unsigned int minPoints, double minTime, double maxTime, double minMomentum,
               bool showElectrons, bool showMuons, bool showGammas, 
               bool showNeutrinos, bool showNeutrons, bool showOthers) :
               TGMainFrame(p, 400, 400), _dataInterface(dataInterface)
  {
    SetCleanup(kDeepCleanup);

    TGHorizontalFrame *subFrame1  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame2  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame3  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame4  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame5  = new TGHorizontalFrame(this,500,20);

    TGLayoutHints *lh = new TGLayoutHints(kLHintsTop,2,1,2,2);
    AddFrame(subFrame1,  lh);
    AddFrame(subFrame2,  lh);
    AddFrame(subFrame3,  lh);
    AddFrame(subFrame4,  lh);

    TGLabel *label1  = new TGLabel(subFrame1,  "Minimum Number of Trajectory Points");
    TGLabel *label2  = new TGLabel(subFrame2,  "Begin Time Window [ns]");
    TGLabel *label3  = new TGLabel(subFrame3,  "End Time Window [ns]");
    TGLabel *label4  = new TGLabel(subFrame4,  "Minumum Momentum [MeV/c]");

    subFrame1->AddFrame(label1, lh);
    subFrame2->AddFrame(label2, lh);
    subFrame3->AddFrame(label3, lh);
    subFrame4->AddFrame(label4, lh);

    _textEntry1 = new TGTextEntry(subFrame1, new TGTextBuffer, 11);
    _textEntry2 = new TGTextEntry(subFrame2, new TGTextBuffer, 12);
    _textEntry3 = new TGTextEntry(subFrame3, new TGTextBuffer, 13);
    _textEntry4 = new TGTextEntry(subFrame4, new TGTextBuffer, 14);

    subFrame1->AddFrame(_textEntry1, lh);
    subFrame2->AddFrame(_textEntry2, lh);
    subFrame3->AddFrame(_textEntry3, lh);
    subFrame4->AddFrame(_textEntry4, lh);

    _textEntry1->Associate(this);
    _textEntry2->Associate(this);
    _textEntry3->Associate(this);
    _textEntry4->Associate(this);

    _textEntry1->SetWidth(100);
    _textEntry2->SetWidth(100);
    _textEntry3->SetWidth(100);
    _textEntry4->SetWidth(100);

    char c[100];
    sprintf(c,"%i",minPoints);   _textEntry1->SetText(c);
    sprintf(c,"%e",minTime);     _textEntry2->SetText(c);
    sprintf(c,"%e",maxTime);     _textEntry3->SetText(c);
    sprintf(c,"%e",minMomentum); _textEntry4->SetText(c);

    _checkButton1 = new TGCheckButton(this,"Show Electrons",21);
    _checkButton2 = new TGCheckButton(this,"Show Muons",22);
    _checkButton3 = new TGCheckButton(this,"Show Gammas",23);
    _checkButton4 = new TGCheckButton(this,"Show Neutrinos",24);
    _checkButton5 = new TGCheckButton(this,"Show Neutrons",25);
    _checkButton6 = new TGCheckButton(this,"Show Other Tracks",26);

    AddFrame(_checkButton1, lh);
    AddFrame(_checkButton2, lh);
    AddFrame(_checkButton3, lh);
    AddFrame(_checkButton4, lh);
    AddFrame(_checkButton5, lh);
    AddFrame(_checkButton6, lh);

    _checkButton1->Associate(this);
    _checkButton2->Associate(this);
    _checkButton3->Associate(this);
    _checkButton4->Associate(this);
    _checkButton5->Associate(this);
    _checkButton6->Associate(this);

    _checkButton1->SetState(showElectrons?kButtonDown:kButtonUp);
    _checkButton2->SetState(showMuons?kButtonDown:kButtonUp);
    _checkButton3->SetState(showGammas?kButtonDown:kButtonUp);
    _checkButton4->SetState(showNeutrinos?kButtonDown:kButtonUp);
    _checkButton5->SetState(showNeutrons?kButtonDown:kButtonUp);
    _checkButton6->SetState(showOthers?kButtonDown:kButtonUp);

    AddFrame(subFrame5,  lh);
    TGTextButton *button1 = new TGTextButton(subFrame5, "&Ok", 101);
    TGTextButton *button2 = new TGTextButton(subFrame5, "&Cancel", 102);
    subFrame5->AddFrame(button1, lh);
    subFrame5->AddFrame(button2, lh);
    button1->Associate(this);
    button2->Associate(this);

    SetWindowName("Filter");
    MapSubwindows();
    MapWindow();
    Layout();

    gClient->WaitFor(this);
  }

  virtual ~FilterDialog()
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
                                            unsigned int minPoints    =atoi(_textEntry1->GetText());
                                            double       minTime      =atof(_textEntry2->GetText());
                                            double       maxTime      =atof(_textEntry3->GetText());
                                            double       minMomentum  =atof(_textEntry4->GetText());
                                            bool showElectrons=(_checkButton1->GetState()==kButtonDown);
                                            bool showMuons    =(_checkButton2->GetState()==kButtonDown);
                                            bool showGammas   =(_checkButton3->GetState()==kButtonDown);
                                            bool showNeutrinos=(_checkButton4->GetState()==kButtonDown);
                                            bool showNeutrons =(_checkButton5->GetState()==kButtonDown);
                                            bool showOthers   =(_checkButton6->GetState()==kButtonDown);
                                            if(_dataInterface) _dataInterface->setTrackFilter(minPoints, minTime, maxTime, minMomentum,
                                                                                              showElectrons, showMuons, showGammas, 
                                                                                              showNeutrinos, showNeutrons, showOthers);
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
#endif /* EventDisplay_src_FilterDialog_h */

