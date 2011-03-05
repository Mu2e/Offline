#include <TSystem.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TColor.h>
#include <TPad.h>
#include <TAxis3D.h>
#include <TView3D.h>
#include <TCanvas.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGTextBuffer.h>
#include <TGTextEntry.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <TGIcon.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TTimer.h>
#include <TText.h>
#include <TBox.h>
#include <TPolyLine.h>

#include "EventDisplayFrame.h"
#include "dict_classes/EventDisplayPad.h"
#include "VirtualShape.h"
#include "DataInterface.h"
#include "ContentSelector.h"

#include "TClassMenuItem.h"
#include "TClass.h"


#include <iostream>

namespace mu2e_eventdisplay
{

EventDisplayFrame::EventDisplayFrame(const TGWindow* p, UInt_t w, UInt_t h) : 
                   TGMainFrame(p, w, h)
{
  int x,y;
  unsigned int width,height;
  gVirtualX->GetWindowSize(gClient->GetRoot()->GetId(),x,y,width,height);
  MoveResize(20,20,width-30,height-70);

  _timer=new TTimer();
  _timer->SetObject(this);
  _timeCurrent=NAN;
  _clock=NULL;
  _isClosed=false;
  _saveAnim=false;
  for(int i=0; i<30; i++)
  {
    _legendText[i]=NULL;
    _legendBox[i]=NULL;
  }
  for(int i=0; i<6; i++)
  {
    _legendParticleText[i]=NULL;
    _legendParticleLine[i]=NULL;
  }

  //bare pointers needed since ROOT manages the following object
  TGHorizontalFrame *mainFrame = new TGHorizontalFrame(this,800,400);
  _mainCanvas = new TRootEmbeddedCanvas("EventDisplayCanvas",mainFrame,GetWidth()-270,GetHeight()-170);
  TGVerticalFrame *subFrame = new TGVerticalFrame(mainFrame,300,400);

  mainFrame->AddFrame(_mainCanvas, new TGLayoutHints(kLHintsTop));
  mainFrame->AddFrame(subFrame, new TGLayoutHints(kLHintsTop));
  AddFrame(mainFrame, new TGLayoutHints(kLHintsTop, 2,2,2,2));

  TGLayoutHints *lh0 = new TGLayoutHints(kLHintsTop,0,0,0,0);
  TGLayoutHints *lh1 = new TGLayoutHints(kLHintsTop,2,1,2,2);

  TGLabel *hitLabel  = new TGLabel(subFrame, "Tracker Hits");
  TGComboBox *hitBox = new TGComboBox(subFrame,10);
  hitBox->Associate(this);
  hitBox->Resize(250,20);
  subFrame->AddFrame(hitLabel, lh1);
  subFrame->AddFrame(hitBox, lh1);

  TGLabel *caloHitLabel  = new TGLabel(subFrame, "Calo Hits");
  TGComboBox *caloHitBox = new TGComboBox(subFrame,11);
  caloHitBox->Associate(this);
  caloHitBox->Resize(250,20);
  subFrame->AddFrame(caloHitLabel, lh1);
  subFrame->AddFrame(caloHitBox, lh1);

  TGLabel *trackLabel  = new TGLabel(subFrame, "Tracks");
  TGListBox *trackBox = new TGListBox(subFrame,12);
  trackBox->Associate(this);
  trackBox->Resize(250,60);
  trackBox->SetMultipleSelections(true);
  subFrame->AddFrame(trackLabel, lh1);
  subFrame->AddFrame(trackBox, lh1);

  _contentSelector=new ContentSelector(hitBox, caloHitBox, trackBox);

  _unhitButton = new TGCheckButton(subFrame,"Show Unhit Straws",31);
  subFrame->AddFrame(_unhitButton, lh1);
  _unhitButton->Associate(this);

  _unhitCrystalsButton = new TGCheckButton(subFrame,"Show Unhit Crystals",36);
  subFrame->AddFrame(_unhitCrystalsButton, lh1);
  _unhitCrystalsButton->Associate(this);

  _supportStructuresButton = new TGCheckButton(subFrame,"Show Tracker Supports, Calo Vanes, Target",32);
  _supportStructuresButton->SetState(kButtonDown);
  subFrame->AddFrame(_supportStructuresButton, lh1);
  _supportStructuresButton->Associate(this);

  _otherStructuresButton = new TGCheckButton(subFrame,"Show Toy DS, CR Steel Shield",37);
  _otherStructuresButton->SetState(kButtonDown);
  subFrame->AddFrame(_otherStructuresButton, lh1);
  _otherStructuresButton->Associate(this);

  _outsideTracksButton = new TGCheckButton(subFrame,"Adjust View to show all Tracks",33);
  subFrame->AddFrame(_outsideTracksButton, lh1);
  _outsideTracksButton->Associate(this);

  _calorimeterViewButton = new TGCheckButton(subFrame,"Adjust View to show Calorimeter",34);
  _calorimeterViewButton->SetState(kButtonDown);
  subFrame->AddFrame(_calorimeterViewButton, lh1);
  _calorimeterViewButton->Associate(this);

  _targetViewButton = new TGCheckButton(subFrame,"Adjust View to show Target",35);
  _targetViewButton->SetState(kButtonDown);
  subFrame->AddFrame(_targetViewButton, lh1);
  _targetViewButton->Associate(this);

  TGHorizontalFrame *subFrameAnim = new TGHorizontalFrame(subFrame,300,15);
  TGTextButton *animButtonStart   = new TGTextButton(subFrameAnim, "Start Animation", 40);
  TGTextButton *animButtonStop    = new TGTextButton(subFrameAnim, "Stop Animation",41);
  TGTextButton *animButtonReset   = new TGTextButton(subFrameAnim, "Reset",42);
  subFrameAnim->AddFrame(animButtonStart, lh1);
  subFrameAnim->AddFrame(animButtonStop, lh1);
  subFrameAnim->AddFrame(animButtonReset, lh1);
  subFrame->AddFrame(subFrameAnim, lh0);
  animButtonStart->Associate(this);
  animButtonStop->Associate(this);
  animButtonReset->Associate(this);

  _repeatAnimationButton = new TGCheckButton(subFrame,"Repeat Animation",42);
  subFrame->AddFrame(_repeatAnimationButton, lh1);
  _repeatAnimationButton->Associate(this);

  TGHorizontalFrame *subFrameSave = new TGHorizontalFrame(subFrame,300,15);
  TGTextButton *saveButton        = new TGTextButton(subFrameSave, "Save", 50);
  TGTextButton *saveAnimButton    = new TGTextButton(subFrameSave, "Save Animation", 51);
  subFrameSave->AddFrame(saveButton, lh1);
  subFrameSave->AddFrame(saveAnimButton, lh1);
  subFrame->AddFrame(subFrameSave, lh0);
  saveButton->Associate(this);
  saveAnimButton->Associate(this);

  _hitColorButton = new TGCheckButton(subFrame,"Use Hit Colors",60);
  _trackColorButton = new TGCheckButton(subFrame,"Use Track Colors",61);
  _backgroundButton = new TGCheckButton(subFrame,"White Background",62);
  subFrame->AddFrame(_hitColorButton, lh1);
  subFrame->AddFrame(_trackColorButton, lh1);
  subFrame->AddFrame(_backgroundButton, lh1);
  _hitColorButton->Associate(this);
  _trackColorButton->Associate(this);
  _backgroundButton->Associate(this);
  _hitColorButton->SetState(kButtonDown);
  _trackColorButton->SetState(kButtonDown);
  _backgroundButton->SetState(kButtonUp);

  _eventInfo = new TGLabel*[4];
  for(int i=0; i<4; i++)
  {
    _eventInfo[i] = new TGLabel(subFrame, "Place Holder for Event Info");
    _eventInfo[i]->SetTextJustify(kTextLeft);
    subFrame->AddFrame(_eventInfo[i], new TGLayoutHints(kLHintsLeft,2,0,2,1));
  }

  TGHorizontalFrame *footLine   = new TGHorizontalFrame(this,800,100);
  _infoCanvas = new TRootEmbeddedCanvas("InfoCanvas",footLine,GetWidth()-600,GetHeight()-450);
  footLine->AddFrame(_infoCanvas, new TGLayoutHints(kLHintsTop));

  TGGroupFrame *zoomangleFrame  = new TGGroupFrame(footLine,"Zoom & Angle");
  TGHorizontalFrame *zoomFrame1 = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGHorizontalFrame *zoomFrame2 = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGHorizontalFrame *angleFrame = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGHorizontalFrame *perspectiveFrame = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGLabel *zoomLabel1  = new TGLabel(zoomFrame1, "minx");
  TGLabel *zoomLabel2  = new TGLabel(zoomFrame1, "mm  miny");
  TGLabel *zoomLabel3  = new TGLabel(zoomFrame1, "mm  minz");
  TGLabel *zoomLabel4  = new TGLabel(zoomFrame1, "mm");
  TGLabel *zoomLabel5  = new TGLabel(zoomFrame2, "maxx");
  TGLabel *zoomLabel6  = new TGLabel(zoomFrame2, "mm  maxy");
  TGLabel *zoomLabel7  = new TGLabel(zoomFrame2, "mm  maxz");
  TGLabel *zoomLabel8  = new TGLabel(zoomFrame2, "mm");
  TGLabel *angleLabel1 = new TGLabel(angleFrame, "phi");
  TGLabel *angleLabel2 = new TGLabel(angleFrame, "°  theta");
  TGLabel *angleLabel3 = new TGLabel(angleFrame, "°  psi");
  TGLabel *angleLabel4 = new TGLabel(angleFrame, "°");
  _minXField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1501);
  _minYField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1502);
  _minZField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1503);
  _maxXField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1504);
  _maxYField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1505);
  _maxZField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1506);
  _phiField   = new TGTextEntry(angleFrame, new TGTextBuffer, 1601);
  _thetaField = new TGTextEntry(angleFrame, new TGTextBuffer, 1602);
  _psiField   = new TGTextEntry(angleFrame, new TGTextBuffer, 1603);
  _perspectiveButton = new TGRadioButton(perspectiveFrame, "perspective", 1700);
  _parallelButton    = new TGRadioButton(perspectiveFrame, "parallel", 1701);
  _minXField->SetWidth(50); 
  _minYField->SetWidth(50); 
  _minZField->SetWidth(50); 
  _maxXField->SetWidth(50); 
  _maxYField->SetWidth(50); 
  _maxZField->SetWidth(50); 
  _phiField->SetWidth(50); 
  _thetaField->SetWidth(50); 
  _psiField->SetWidth(50); 
  TGTextButton *setRangeButton = new TGTextButton(zoomangleFrame, "Set &Range", 1500);
  TGTextButton *setAngleButton = new TGTextButton(zoomangleFrame, "Set &Angle", 1600);
  zoomFrame1->AddFrame(zoomLabel1, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(_minXField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(zoomLabel2, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(_minYField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(zoomLabel3, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(_minZField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame1->AddFrame(zoomLabel4, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(zoomLabel5, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(_maxXField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(zoomLabel6, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(_maxYField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(zoomLabel7, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(_maxZField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame2->AddFrame(zoomLabel8, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(angleLabel1, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(_phiField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(angleLabel2, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(_thetaField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(angleLabel3, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(_psiField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  angleFrame->AddFrame(angleLabel4, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  perspectiveFrame->AddFrame(_perspectiveButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  perspectiveFrame->AddFrame(_parallelButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(zoomFrame1, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(zoomFrame2, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(setRangeButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(angleFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(setAngleButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(perspectiveFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  footLine->AddFrame(zoomangleFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  _minXField->Associate(this);
  _minYField->Associate(this);
  _minZField->Associate(this);
  _maxXField->Associate(this);
  _maxYField->Associate(this);
  _maxZField->Associate(this);
  _phiField->Associate(this);
  _thetaField->Associate(this);
  _psiField->Associate(this);
  setRangeButton->Associate(this);
  setAngleButton->Associate(this);
  _perspectiveButton->Associate(this);
  _parallelButton->Associate(this);
  _perspectiveButton->SetState(kButtonDown);
  _parallelButton->SetState(kButtonUp);

  TGVerticalFrame *innerFrame1   = new TGVerticalFrame(footLine,100,400);

  TGGroupFrame *optionsFrame     = new TGGroupFrame(innerFrame1,"Options");
  TGHorizontalFrame *filterFrame = new TGHorizontalFrame(optionsFrame,500,50);
  TGHorizontalFrame *jumpFrame   = new TGHorizontalFrame(optionsFrame,500,50);
  TGLabel *minHitLabel      = new TGLabel(filterFrame, "minimum hits");
  TGLabel *eventToFindLabel = new TGLabel(jumpFrame, "jump to event number");
  _minHitField = new TGTextEntry(filterFrame, new TGTextBuffer, 1101);
  _eventToFindField = new TGTextEntry(jumpFrame, new TGTextBuffer, 1103);
  _minHits=0;  
  _minHitField->SetWidth(50); 
  _minHitField->SetText("0");
  _findEvent=false;
  _eventToFind=0;
  _eventToFindField->SetWidth(50); 
  _eventToFindField->SetText("");
  TGTextButton *applyButton = new TGTextButton(filterFrame, "&Apply", 1100);
  TGTextButton *goButton    = new TGTextButton(jumpFrame, "&Go", 1102);
  filterFrame->AddFrame(minHitLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  filterFrame->AddFrame(_minHitField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  filterFrame->AddFrame(applyButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  jumpFrame->AddFrame(eventToFindLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  jumpFrame->AddFrame(_eventToFindField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  jumpFrame->AddFrame(goButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,3,0,3,0));
  optionsFrame->AddFrame(filterFrame, new TGLayoutHints(kLHintsLeft,3,0,0,0));
  optionsFrame->AddFrame(jumpFrame, new TGLayoutHints(kLHintsLeft,3,0,0,0));
  innerFrame1->AddFrame(optionsFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  TGHorizontalFrame *navigationFrame = new TGHorizontalFrame(innerFrame1,100,200);
  TGTextButton *quitButton         = new TGTextButton(navigationFrame, "&Quit", 1001);
  TGTextButton *nextButton         = new TGTextButton(navigationFrame, "&Next", 1111);
  navigationFrame->AddFrame(quitButton, new TGLayoutHints(kLHintsLeft,10,0,10,0));
  navigationFrame->AddFrame(nextButton, new TGLayoutHints(kLHintsLeft,10,0,10,0));
  innerFrame1->AddFrame(navigationFrame, new TGLayoutHints(kLHintsLeft,10,0,10,0));

  quitButton->Associate(this);
  nextButton->Associate(this);
  _minHitField->Associate(this);
  _eventToFindField->Associate(this);
  applyButton->Associate(this);
  goButton->Associate(this);

  std::string logoFileName=getenv("MU2E_BASE_RELEASE");
  logoFileName.append("/EventDisplay/src/logo_small.png");
  const TGPicture *logo = gClient->GetPicture(logoFileName.c_str());
  TGIcon *icon = new TGIcon(navigationFrame, logo, 50, 50);
  navigationFrame->AddFrame(icon, new TGLayoutHints(kLHintsLeft,20,0,0,0));

  footLine->AddFrame(innerFrame1, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  AddFrame(footLine, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  MapSubwindows();
  SetWindowName("Mu2e Event Display");
  MapWindow();

  _mainCanvas->GetCanvas()->cd();
  _mainPad = new EventDisplayPad("mainPad","Detector", 0, 0, 1, 1, 5,1,1);  
  _mainPad->setEventDisplayFrame(this);
  _mainPad->SetFillColor(1);
  _mainPad->Draw();

  _infoCanvas->GetCanvas()->cd();
  _infoPad = new TPad("infoPad","InfoField", 0, 0, 1, 1, 5,1,1);  
  _infoPad->SetFillColor(0);
  _infoPad->Draw();

  for(int i=0; i<20; i++)
  {
    float r,g,b;
    TColor *c;
    TColor::HLS2RGB(i*360/20,.5,.5,r,g,b);
    if(!gROOT->GetColor(i+2000)) c = new TColor(i+2000,r,g,b);
  }

  _mainPad->cd();
  _dataInterface = boost::shared_ptr<DataInterface>(new DataInterface(this));
}

EventDisplayFrame::~EventDisplayFrame()
{
  // TODO
  // if(timer) delete timer;
  // Cleanup();
}

Bool_t EventDisplayFrame::HandleConfigureNotify(Event_t *event)
{
// This is a modified version of the function from TGFrame.cxx
   if ((event->fWidth != fWidth) || (event->fHeight != fHeight)) 
   {
      fWidth  = event->fWidth;
      fHeight = event->fHeight;
      _mainCanvas->SetWidth(fWidth-270);
      _mainCanvas->SetHeight(fHeight-170);
      _infoCanvas->SetWidth(fWidth-600);
      _infoCanvas->SetHeight(fHeight-450);
      Layout();
   }
   return kTRUE;
}

void EventDisplayFrame::fillZoomAngleFields()
{
  if(_mainPad->GetView()==NULL) return;
  char c[100];
  double min[3], max[3];
  _mainPad->GetView()->GetRange(min,max);
  sprintf(c,"%.0f",min[0]); _minXField->SetText(c);
  sprintf(c,"%.0f",min[1]); _minYField->SetText(c);
  sprintf(c,"%.0f",min[2]); _minZField->SetText(c);
  sprintf(c,"%.0f",max[0]); _maxXField->SetText(c);
  sprintf(c,"%.0f",max[1]); _maxYField->SetText(c);
  sprintf(c,"%.0f",max[2]); _maxZField->SetText(c);
  sprintf(c,"%.0f",_mainPad->GetView()->GetLongitude()); _phiField->SetText(c);
  sprintf(c,"%.0f",_mainPad->GetView()->GetLatitude()); _thetaField->SetText(c);
  sprintf(c,"%.0f",_mainPad->GetView()->GetPsi()); _psiField->SetText(c);
  if(_mainPad->GetView()->IsPerspective())
  {
    _perspectiveButton->SetState(kButtonDown);
    _parallelButton->SetState(kButtonUp);
  }
  else
  {
    _perspectiveButton->SetState(kButtonUp);
    _parallelButton->SetState(kButtonDown);
  }
}

void EventDisplayFrame::fillGeometry()
{
  _mainPad->cd();
  _dataInterface->fillGeometry();
  DataInterface::spaceminmax m=_dataInterface->getSpaceBoundary(false, true, false);
  _mainPad->GetView()->SetRange(m.minx,m.miny,m.minz,m.maxx,m.maxy,m.maxz);
  _mainPad->GetView()->AdjustScales();
  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::setEvent(const edm::Event& event, bool firstLoop)
{
  char eventInfoText[50];
  sprintf(eventInfoText,"Event #: %i",event.id().event());
  _eventInfo[0]->SetText(eventInfoText);
  sprintf(eventInfoText,"Run #: %i",event.id().run());
  _eventInfo[1]->SetText(eventInfoText);
  this->Layout();

  _contentSelector->setAvailableCollections(event);
  if(firstLoop) _contentSelector->firstLoop();
  fillEvent(firstLoop);

  gApplication->Run(true);
}

void EventDisplayFrame::fillEvent(bool firstLoop)
{
  _findEvent=false;
  _mainPad->cd();
  _dataInterface->fillEvent(_contentSelector);
  _dataInterface->useHitColors(_hitColorButton->GetState()==kButtonDown,
                               _backgroundButton->GetState()==kButtonDown);
  _dataInterface->useTrackColors(_trackColorButton->GetState()==kButtonDown,
                                 _backgroundButton->GetState()==kButtonDown);
  updateHitLegend(_hitColorButton->GetState()==kButtonDown);
  updateTrackLegend(_trackColorButton->GetState()==kButtonDown);

  //set zoom only if "all tracks" option is on, or if it is the first event
  if(_outsideTracksButton->GetState()==kButtonDown || firstLoop) 
  {
    DataInterface::spaceminmax m=_dataInterface->getSpaceBoundary(_targetViewButton->GetState()==kButtonDown,
                                                   _calorimeterViewButton->GetState()==kButtonDown,
                                                   _outsideTracksButton->GetState()==kButtonDown);
    _mainPad->GetView()->SetRange(m.minx,m.miny,m.minz,m.maxx,m.maxy,m.maxz);
    _mainPad->GetView()->AdjustScales();
  }

  char eventInfoText[50];
  sprintf(eventInfoText,"Number of hit straws: %i",_dataInterface->getNumberHits());
  _eventInfo[2]->SetText(eventInfoText);
  sprintf(eventInfoText,"Number of hit calorimeter crystals: %i",_dataInterface->getNumberCrystalHits());
  _eventInfo[3]->SetText(eventInfoText);
  this->Layout();

  drawEverything();
}

void EventDisplayFrame::updateHitLegend(bool draw)
{
  for(int i=0; i<21; i++)
  {
    if(_legendBox[i]!=NULL) delete _legendBox[i];
    if(_legendText[i]!=NULL) delete _legendText[i];
    _legendBox[i]=NULL;
    _legendText[i]=NULL;
  }

  if(draw)
  {
    for(int i=0; i<21; i++)
    {
      if(i<20)
      {
        _legendBox[i]=new TBox(0.6,0.557+i*0.02,0.7,0.577+i*0.02);
        _legendBox[i]->SetFillColor(i+2000);
        _legendBox[i]->Draw();
      }
      if(i%4==0)
      {
        double mint=_dataInterface->getHitsTimeBoundary().mint;
        double maxt=_dataInterface->getHitsTimeBoundary().maxt;
        double t=i*(maxt-mint)/20.0+mint;
        char s[50];
        sprintf(s,"%+.3e ns",t);
        _legendText[i]=new TText(0.72,0.54+i*0.02,s);
        _legendText[i]->SetTextColor(kGray);
        _legendText[i]->SetTextSize(0.025);
        _legendText[i]->Draw();
      }
    }
  }
}

void EventDisplayFrame::updateTrackLegend(bool draw)
{
  for(int i=0; i<6; i++)
  {
    if(_legendParticleLine[i]!=NULL) delete _legendParticleLine[i];
    if(_legendParticleText[i]!=NULL) delete _legendParticleText[i];
    _legendParticleLine[i]=NULL;
    _legendParticleText[i]=NULL;
  }

  if(draw)
  {
    for(int i=0; i<6; i++)
    {
      _legendParticleLine[i]=new TPolyLine();
      _legendParticleLine[i]->SetPoint(0, 0.6,0.45-i*0.05);
      _legendParticleLine[i]->SetPoint(1, 0.7,0.45-i*0.05);
      _legendParticleLine[i]->Draw();
      _legendParticleText[i]=new TText(0.72,0.44-i*0.05,"");
      _legendParticleText[i]->SetTextColor(kGray);
      _legendParticleText[i]->SetTextSize(0.025);
      _legendParticleText[i]->Draw();
    }
    _legendParticleLine[0]->SetLineColor(2);
    _legendParticleLine[1]->SetLineColor(3);
    _legendParticleLine[2]->SetLineColor(4); 
    _legendParticleLine[3]->SetLineColor(6);
    _legendParticleLine[4]->SetLineColor(28);
    _legendParticleLine[5]->SetLineColor(kGray);
    _legendParticleText[0]->SetTitle("e+, e-");
    _legendParticleText[1]->SetTitle("mu+, mu-");
    _legendParticleText[2]->SetTitle("gamma");
    _legendParticleText[3]->SetTitle("n0");
    _legendParticleText[4]->SetTitle("neutrinos");
    _legendParticleText[5]->SetTitle("other particles");
  }
}

bool EventDisplayFrame::isClosed() const
{
  return _isClosed;
}

bool EventDisplayFrame::getSelectedHitsName(std::string &className, 
                                            std::string &moduleLabel, 
                                            std::string &productInstanceName) const
{
  return _contentSelector->getSelectedHitsName(className, moduleLabel, productInstanceName);
}

int EventDisplayFrame::getMinimumHits() const
{
  return _minHits;
}

int EventDisplayFrame::getEventToFind(bool &findEvent) const
{
  findEvent=_findEvent;
  return _eventToFind;
}

void EventDisplayFrame::CloseWindow()
{
  _isClosed=true;
  _timer->Stop();
  _timeCurrent=NAN;
  TGMainFrame::CloseWindow();
  gApplication->Terminate();
}

Bool_t EventDisplayFrame::ProcessMessage(Long_t msg, Long_t param1, Long_t param2)
{
  switch (GET_MSG(msg))
  {
    case kC_COMMAND: 
      switch (GET_SUBMSG(msg))
      {
        case kCM_BUTTON: if(param1==1001) CloseWindow();
                         if(param1==40) prepareAnimation();
                         if(param1==41) 
                         {
                           _timer->Stop(); 
                           _timeCurrent=NAN;
                           if(_saveAnim) combineAnimFiles();
                         }
                         if(param1==42) 
                         {
                           _timer->Stop(); 
                           _timeCurrent=NAN;
                           drawEverything();
                         }
                         if(param1==1111)
                         {
                           _timer->Stop();
                           _timeCurrent=NAN;
                           _contentSelector->setSelectedHitsName();
                           gApplication->Terminate();
                         }
                         if(param1==1100)
                         {
                           _minHits=atoi(_minHitField->GetText());
                         }
                         if(param1==1102)
                         {
                           _eventToFind=atoi(_eventToFindField->GetText()); 
                           _findEvent=true;
                           _timer->Stop();
                           _timeCurrent=NAN;
                           gApplication->Terminate();
                         }
                         if(param1==50 || param1==51)
                         {
                           const char *fileType[]={"Gif files","*.gif",0,0};
                           TGFileInfo fileInfo;
                           fileInfo.fFileTypes = fileType;
                           TGFileDialog *fileDialog;
                           fileDialog=new TGFileDialog(gClient->GetRoot(), gClient->GetRoot(), kFDSave, &fileInfo); //ROOT takes care of deleting this
                           if(!fileInfo.fFilename) break;
                           char f[strlen(fileInfo.fFilename)+5];
                           strcpy(f,fileInfo.fFilename);
                           if(strcmp(f+strlen(f)-4,".gif")!=0) strcat(f,".gif");
                           if(param1==50) _mainPad->SaveAs(f);
                           if(param1==51)
                           {
                             _saveAnim=true;
                             _saveAnimFile.assign(f,strlen(f)-4);
                             prepareAnimation();
                           }
                         }
                         if(param1==1500)
                         {
                           double min[3],max[3];
                           min[0]=atof(_minXField->GetText());
                           min[1]=atof(_minYField->GetText());
                           min[2]=atof(_minZField->GetText());
                           max[0]=atof(_maxXField->GetText());
                           max[1]=atof(_maxYField->GetText());
                           max[2]=atof(_maxZField->GetText());
                           if(min[0]<max[0] && min[1]<max[1] && min[2]<max[2])
                             _mainPad->GetView()->SetRange(min,max);
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         if(param1==1600)
                         {
                           double phi=atof(_phiField->GetText());
                           double theta=atof(_thetaField->GetText());
                           double psi=atof(_psiField->GetText());
                           int irep=0;
                           _mainPad->GetView()->SetView(phi,theta,psi,irep);
                           _mainPad->SetPhi(-90-phi);
                           _mainPad->SetTheta(90-theta);
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         break;
   case kCM_RADIOBUTTON: if(param1==1700)
                         {
                           _parallelButton->SetState(kButtonUp);
                           _perspectiveButton->SetState(kButtonDown);
                           _mainPad->GetView()->SetPerspective();
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         if(param1==1701)
                         {
                           _perspectiveButton->SetState(kButtonUp);
                           _parallelButton->SetState(kButtonDown);
                           _mainPad->GetView()->SetParallel();
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         break;
   case kCM_CHECKBUTTON: if(param1==31)
                         {
                           _mainPad->cd();
                           if(_unhitButton->GetState()==kButtonDown)
                           {
                             _dataInterface->makeStrawsVisibleBeforeStart(true);
                           }
                           else
                           {
                             _dataInterface->makeStrawsVisibleBeforeStart(false);
                           }
                           drawEverything();
                         }
                         if(param1==32)
                         {
                           _mainPad->cd();
                           if(_supportStructuresButton->GetState()==kButtonDown)
                           {
                             _dataInterface->makeSupportStructuresVisible(true);
                           }
                           else
                           {
                             _dataInterface->makeSupportStructuresVisible(false);
                           }
                           drawEverything();
                         }
                         if(param1==37)
                         {
                           _mainPad->cd();
                           if(_otherStructuresButton->GetState()==kButtonDown)
                           {
                             _dataInterface->makeOtherStructuresVisible(true);
                           }
                           else
                           {
                             _dataInterface->makeOtherStructuresVisible(false);
                           }
                           drawEverything();
                         }
                         if(param1==36)
                         {
                           _mainPad->cd();
                           if(_unhitCrystalsButton->GetState()==kButtonDown)
                           {
                             _dataInterface->makeCrystalsVisibleBeforeStart(true);
                           }
                           else
                           {
                             _dataInterface->makeCrystalsVisibleBeforeStart(false);
                           }
                           drawEverything();
                         }
                         if(param1==33 || param1==34 || param1==35)
                         {
                           _mainPad->cd();
                           DataInterface::spaceminmax m=_dataInterface->getSpaceBoundary(
                                                 _targetViewButton->GetState()==kButtonDown,
                                                 _calorimeterViewButton->GetState()==kButtonDown,
                                                 _outsideTracksButton->GetState()==kButtonDown);
                           _mainPad->GetView()->SetRange(m.minx,m.miny,m.minz,m.maxx,m.maxy,m.maxz);
                           _mainPad->GetView()->AdjustScales();
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         if(param1==60)
                         {
                           _mainPad->cd();
                           _dataInterface->useHitColors(_hitColorButton->GetState()==kButtonDown,
                                                        _backgroundButton->GetState()==kButtonDown);
                           updateHitLegend(_hitColorButton->GetState()==kButtonDown);
                           if(isnan(_timeCurrent)) drawEverything();
                           else drawSituation();
                         }
                         if(param1==61)
                         {
                           _mainPad->cd();
                           _dataInterface->useTrackColors(_trackColorButton->GetState()==kButtonDown,
                                                          _backgroundButton->GetState()==kButtonDown);
                           updateTrackLegend(_trackColorButton->GetState()==kButtonDown);
                           if(isnan(_timeCurrent)) drawEverything();
                           else drawSituation();
                         }
                         if(param1==62)
                         {
                           _mainPad->cd();
                           if(_backgroundButton->GetState()==kButtonDown) _mainPad->SetFillColor(0);
                           else _mainPad->SetFillColor(1);
                           _dataInterface->useHitColors(_hitColorButton->GetState()==kButtonDown,
                                                        _backgroundButton->GetState()==kButtonDown);
                           _dataInterface->useTrackColors(_trackColorButton->GetState()==kButtonDown,
                                                          _backgroundButton->GetState()==kButtonDown);
                           if(isnan(_timeCurrent)) drawEverything();
                           else drawSituation();
                         }
                         break;
  case kCM_COMBOBOX : if(param1==10) fillEvent();
                      if(param1==11) fillEvent();
                      break;
  case kCM_LISTBOX : if(param1==12) fillEvent();
                     break;
      }
      break;
  }
  return kTRUE;
}
                                        
void EventDisplayFrame::prepareAnimation()
{
  _timer->Stop();           //needed if an animation is already running
  if(_clock) {delete _clock; _clock=NULL;}
  _mainPad->cd();
  _dataInterface->startComponents();
  _mainPad->Modified();
  _mainPad->Update();
 
  if(_outsideTracksButton->GetState()==kButtonDown)
  {
    _timeStart=_dataInterface->getTracksTimeBoundary().mint; 
    _timeStop=_dataInterface->getTracksTimeBoundary().maxt;
  }
  else
  {
    _timeStart=_dataInterface->getHitsTimeBoundary().mint; 
    _timeStop=_dataInterface->getHitsTimeBoundary().maxt; 
  }
  if(isnan(_timeStart) || isnan(_timeStop)) return;
  double diff=_timeStop-_timeStart;
  _timeStart-=diff*0.05;
  _timeStop+=diff*0.05;
  _timeCurrent=_timeStart;
  _saveAnimCounter=0;
  _timer->Start(100,kFALSE);
}


Bool_t EventDisplayFrame::HandleTimer(TTimer *)
{
  _timer->Reset();
  _timeCurrent+=(_timeStop-_timeStart)/200;
  drawSituation();
  if(_saveAnim) 
  {
    _saveAnimCounter++;
    if(_saveAnimCounter%3==0) //save only every 3rd gif to make final file smaller
    {
      char c[_saveAnimFile.length()+15];
      sprintf(c,"%s_tmp_%04i.gif",_saveAnimFile.c_str(),_saveAnimCounter);
      _mainPad->SaveAs(c);
    }
  }
  if(_timeCurrent>=_timeStop) 
  {
    _timer->Stop();
    _timeCurrent=NAN;
    if(_saveAnim) combineAnimFiles();
    if(_repeatAnimationButton->GetState()==kButtonDown) prepareAnimation();
  }
  return kTRUE;
}

void EventDisplayFrame::combineAnimFiles()
{
  _saveAnim=false;
  char c[2*_saveAnimFile.length()+100];
  sprintf(c,"convert -delay 50 -loop 0 %s_tmp_*.gif %s.gif",_saveAnimFile.c_str(),_saveAnimFile.c_str());
  gSystem->Exec(c);
  sprintf(c,"rm -f %s_tmp_*.gif",_saveAnimFile.c_str());
  gSystem->Exec(c);
}

void EventDisplayFrame::drawSituation()
{
  _mainPad->cd();
  _dataInterface->updateComponents(_timeCurrent);
  if(!_clock)
  {
    char timeText[50];
    sprintf(timeText,"%+.4e ns",_timeCurrent);
    _clock = new TText(0.52,-0.9,timeText);
    _clock->SetTextColor(5);
    _clock->Draw("same");
  }
  else
  {
    char timeText[50];
    sprintf(timeText,"%+.4e ns",_timeCurrent);
    _clock->SetTitle(timeText);
  }

  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::drawEverything()
{
  _mainPad->cd();
  if(_clock) {delete _clock; _clock=NULL;}
  if(TAxis3D::GetPadAxis(_mainPad)==NULL) _mainPad->GetView()->ShowAxis();
  TAxis3D::GetPadAxis(_mainPad)->SetLabelSize(0.025); 
  _mainPad->Modified();
  _mainPad->Update();
  _dataInterface->updateComponents(NAN);
  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::showInfo(TObject *o)  //ROOT accepts only bare pointers here
{
  _infoPad->cd();  
  _infoPad->Clear();
  (dynamic_cast<ComponentInfo*>(o))->showInfo();
  _infoPad->Modified();
  _infoPad->Update();
  _mainPad->cd();  
}

}
