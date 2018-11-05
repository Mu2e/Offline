#ifndef EventDisplay_src_FilterDialog_h
#define EventDisplay_src_FilterDialog_h

#include <TGFrame.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGTextEntry.h>
#include "EventDisplay/src/DataInterface.h"
#include "EventDisplay/src/ContentSelector.h"

namespace mu2e_eventdisplay
{

class FilterDialog : public TGMainFrame
{
  boost::shared_ptr<DataInterface> _dataInterface;
  boost::shared_ptr<ContentSelector> _contentSelector;

  TGTextEntry   *_textEntry1;
  TGTextEntry   *_textEntry2;
  TGTextEntry   *_textEntry3;
  TGTextEntry   *_textEntry4;
  TGCheckButton *_checkButton[6];
  TGCheckButton *_checkButtonFlags[32];
  TGComboBox    *_hitFlagBox;

  FilterDialog();
  FilterDialog(const FilterDialog &);
  FilterDialog& operator=(const FilterDialog &);

  public:
  FilterDialog(const TGWindow* p, 
               boost::shared_ptr<DataInterface> dataInterface, 
               boost::shared_ptr<ContentSelector> contentSelector) : 
               TGMainFrame(p, 500, 700), _dataInterface(dataInterface), _contentSelector(contentSelector)
  {
    SetCleanup(kDeepCleanup);

    unsigned int minPoints=0;
    double       minTime=NAN;
    double       maxTime=NAN;
    double       minMomentum=0;
    bool showElectrons=true;
    bool showMuons=true;
    bool showGammas=true;
    bool showNeutrinos=true;
    bool showNeutrons=true;
    bool showOthers=true;
    
    mu2e::StrawHitFlag hitFlagSelection;

    if(_dataInterface) _dataInterface->getFilterValues(minPoints, minTime, maxTime, minMomentum,
                                                       showElectrons, showMuons, showGammas, 
                                                       showNeutrinos, showNeutrons, showOthers, 
                                                       hitFlagSelection);

    TGHorizontalFrame *subFrame1  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame2  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame3  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame4  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame5  = new TGHorizontalFrame(this,500,20);
    TGHorizontalFrame *subFrame6  = new TGHorizontalFrame(this,500,20);

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

    _checkButton[0] = new TGCheckButton(this,"Show Electrons",21);
    _checkButton[1] = new TGCheckButton(this,"Show Muons",22);
    _checkButton[2] = new TGCheckButton(this,"Show Gammas",23);
    _checkButton[3] = new TGCheckButton(this,"Show Neutrinos",24);
    _checkButton[4] = new TGCheckButton(this,"Show Neutrons",25);
    _checkButton[5] = new TGCheckButton(this,"Show Other Tracks",26);

    for(int i=0; i<6; i++)
    {
      AddFrame(_checkButton[i], lh);
      _checkButton[i]->Associate(this);
    }

    _checkButton[0]->SetState(showElectrons?kButtonDown:kButtonUp);
    _checkButton[1]->SetState(showMuons?kButtonDown:kButtonUp);
    _checkButton[2]->SetState(showGammas?kButtonDown:kButtonUp);
    _checkButton[3]->SetState(showNeutrinos?kButtonDown:kButtonUp);
    _checkButton[4]->SetState(showNeutrons?kButtonDown:kButtonUp);
    _checkButton[5]->SetState(showOthers?kButtonDown:kButtonUp);

    AddFrame(subFrame5, lh);
    TGGroupFrame *hitFlagGroupFrame  = new TGGroupFrame(subFrame5,"Hit Flags");
    subFrame5->AddFrame(hitFlagGroupFrame); 
 
    _hitFlagBox = new TGComboBox(hitFlagGroupFrame,30);
    _hitFlagBox->Resize(250,20);
    _hitFlagBox->Associate(this);
    hitFlagGroupFrame->AddFrame(_hitFlagBox, lh);

    if(_contentSelector)
    {
      const std::vector<ContentSelector::entryStruct> &hitFlagEntries = _contentSelector->getStrawHitFlagEntries();
      const std::string &selectedEntry = _contentSelector->getSelectedStrawHitFlagEntry();
      for(unsigned int i=0; i<hitFlagEntries.size(); i++)
      {
        _hitFlagBox->AddEntry(hitFlagEntries[i].entryText.c_str(), hitFlagEntries[i].entryID);
        if(hitFlagEntries[i].entryText.compare(selectedEntry)==0) _hitFlagBox->Select(hitFlagEntries[i].entryID);
      }
      _hitFlagBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
    }


    TGHorizontalFrame *hitFlagFrame0 = new TGHorizontalFrame(hitFlagGroupFrame,600,500);
    TGVerticalFrame *hitFlagFrame1 = new TGVerticalFrame(hitFlagFrame0,300,500);
    TGVerticalFrame *hitFlagFrame2 = new TGVerticalFrame(hitFlagFrame0,300,500);
    hitFlagGroupFrame->AddFrame(hitFlagFrame0); 
    hitFlagFrame0->AddFrame(hitFlagFrame1); 
    hitFlagFrame0->AddFrame(hitFlagFrame2);

    std::string flagnames[16]={"Stereo","EnergySel","RadSel","TimeSel","","",
                               "Delta","Isolated","Outlier","Other","","",
                               "CaloSel","","",""};
    for(int i=0; i<16; i++)
    {
      _checkButtonFlags[i] = new TGCheckButton(hitFlagFrame1,flagnames[i].c_str(),1000+i);
      hitFlagFrame1->AddFrame(_checkButtonFlags[i], lh);
    }
    for(int i=16; i<32; i++)
    {
      char trackname[10];
      sprintf(trackname,"Track%i",i-16);
      _checkButtonFlags[i] = new TGCheckButton(hitFlagFrame2,trackname,1000+i);
      hitFlagFrame2->AddFrame(_checkButtonFlags[i], lh);
    }
    for(int i=0; i<32; i++)
    {
      _checkButtonFlags[i]->Associate(this);
      mu2e::StrawHitFlagDetail::bit_type b=static_cast<mu2e::StrawHitFlagDetail::bit_type>(i);
      _checkButtonFlags[i]->SetState(hitFlagSelection.hasAnyProperty(b)?kButtonDown:kButtonUp);
    }

    AddFrame(subFrame6,  lh);
    TGTextButton *button1 = new TGTextButton(subFrame6, "&Ok", 101);
    TGTextButton *button2 = new TGTextButton(subFrame6, "&Cancel", 102);
    subFrame6->AddFrame(button1, lh);
    subFrame6->AddFrame(button2, lh);
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
                                            bool showElectrons=(_checkButton[0]->GetState()==kButtonDown);
                                            bool showMuons    =(_checkButton[1]->GetState()==kButtonDown);
                                            bool showGammas   =(_checkButton[2]->GetState()==kButtonDown);
                                            bool showNeutrinos=(_checkButton[3]->GetState()==kButtonDown);
                                            bool showNeutrons =(_checkButton[4]->GetState()==kButtonDown);
                                            bool showOthers   =(_checkButton[5]->GetState()==kButtonDown);
                                            mu2e::StrawHitFlag hitFlagSetting;
                                            for(int j=0; j<32; j++)
                                            {
                                              mu2e::StrawHitFlagDetail::bit_type b=static_cast<mu2e::StrawHitFlagDetail::bit_type>(j);
                                              if(_checkButtonFlags[j]->GetState()==kButtonDown) hitFlagSetting.merge(b);
                                            }
                                            if(_dataInterface) 
                                            {
                                              _dataInterface->setFilterValues(minPoints, minTime, maxTime, minMomentum,
                                                                              showElectrons, showMuons, showGammas, 
                                                                              showNeutrinos, showNeutrons, showOthers,
                                                                              hitFlagSetting);
                                              TGTextLBEntry *selectedEntry=dynamic_cast<TGTextLBEntry*>(_hitFlagBox->GetSelectedEntry());
                                              std::string selectedEntryText = selectedEntry->GetText()->GetString();
                                              if(selectedEntry && _contentSelector)
                                              {
                                                _contentSelector->setSelectedStrawHitFlagEntry(selectedEntryText);
                                              }
                                            }
                                            CloseWindow();
                                          }
                                          break;
                       }
                       break;
    }
    return kTRUE;
  }

  ClassDef(FilterDialog,1);
};

}
#endif /* EventDisplay_src_FilterDialog_h */

