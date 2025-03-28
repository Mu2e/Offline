#include <TApplication.h>
#include <TAxis3D.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGIcon.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGTextBuffer.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TPad.h>
#include <TPolyLine.h>
#include <TROOT.h>
#include <TRootEmbeddedCanvas.h>
#include <TSystem.h>
#include <TText.h>
#include <TTimer.h>
#include <TView3D.h>
using namespace std;
#include "Offline/EventDisplay/src/TrackColorSelector.h"
#include "Offline/EventDisplay/src/ContentSelector.h"
#include "Offline/EventDisplay/src/SaveDialogManager.h"
#include "Offline/EventDisplay/src/RootFileManager.h"
#include "Offline/EventDisplay/src/DataInterface.h"
#include "Offline/EventDisplay/src/EventDisplayFrame.h"
#include "Offline/EventDisplay/src/FilterDialog.h"
#include "Offline/EventDisplay/src/SetupDialog.h"
#include "Offline/EventDisplay/src/VirtualShape.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayViewSetup.h"
#include "Offline/EventDisplay/src/dict_classes/ComponentInfoContainer.h"
#include "Offline/EventDisplay/src/dict_classes/HistDraw.h"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"

#include "fhiclcpp/ParameterSet.h"

#include "TGGC.h"
#include "TGFont.h"
#include "TClass.h"
#include "TClassMenuItem.h"

#include <iostream>

namespace mu2e_eventdisplay
{

EventDisplayFrame::EventDisplayFrame(const TGWindow* p, UInt_t w, UInt_t h, fhicl::ParameterSet const &pset) :
  TGMainFrame(p, w, h),
  _wideband(pset.get<bool>("wideband",false)),
  _extracted(pset.get<bool>("extracted",false)),
  _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
  _physicalVolumesMultiLabel(pset.get<std::string>("physicalVolumesMultiLabel","compressPV")),
  _protonBunchTimeLabel(pset.get<std::string>("protonBunchTimeTag","EWMProducer")),
  _kalStepSize(pset.get<double>("kalSeedStepSize"))
{
  SetCleanup(kDeepCleanup);
  Move(20,20);

  FontStruct_t buttonfont = gClient->GetFontByName("-*-helvetica-medium-r-*-*-8-*-*-*-*-*-iso8859-1");
  GCValues_t gval;
  gval.fMask = kGCForeground | kGCFont;
  gval.fFont = gVirtualX->GetFontHandle(buttonfont);
  gClient->GetColorByName("black", gval.fForeground);
  GContext_t buttoncontext = gVirtualX->CreateGC(gClient->GetRoot()->GetId(), &gval);

  _timer=new TTimer();
  _timer->SetObject(this);
  _timeCurrent=NAN;
  _eventNumberText=nullptr;
  _subrunNumberText=nullptr;
  _runNumberText=nullptr;
  _clock=nullptr;
  _isClosed=false;
  _saveAnim=false;
  for(int i=0; i<30; i++)
  {
    _legendText[i]=nullptr;
    _legendBox[i]=nullptr;
  }
  for(int i=0; i<30; i++)
  {
    _legendParticleGroup[i]=nullptr;
    _legendParticleText[i]=nullptr;
    _legendParticleLine[i]=nullptr;
  }

  initSetup();

  //bare pointers needed since ROOT manages the following object
  _mainFrame = new TGHorizontalFrame(this,GetWidth()-270,GetHeight()*0.7);
  _mainCanvas = new TRootEmbeddedCanvas("EventDisplayCanvas",_mainFrame,GetWidth()-270,GetHeight()*0.7);
  _subFrame = new TGVerticalFrame(_mainFrame,270,GetHeight()*0.7);
  _subFrame->SetMaxHeight(GetHeight()*0.7);

  _mainFrame->AddFrame(_mainCanvas, new TGLayoutHints(kLHintsTop));
  _mainFrame->AddFrame(_subFrame, new TGLayoutHints(kLHintsTop));
  AddFrame(_mainFrame, new TGLayoutHints(kLHintsTop, 2,2,2,2));

  TGLayoutHints *lh0 = new TGLayoutHints(kLHintsTop,0,0,0,0);
  TGLayoutHints *lh1 = new TGLayoutHints(kLHintsTop,2,1,2,2);

  TGLabel *hitLabel  = new TGLabel(_subFrame, "Tracker Hits");
  TGComboBox *hitBox = new TGComboBox(_subFrame,10);
  hitBox->Associate(this);
  hitBox->Resize(250,20);
  _subFrame->AddFrame(hitLabel, lh1);
  _subFrame->AddFrame(hitBox, lh1);

  TGLabel *caloHitLabel  = new TGLabel(_subFrame, "Calo Hits");
  TGComboBox *caloHitBox = new TGComboBox(_subFrame,11);
  caloHitBox->Associate(this);
  caloHitBox->Resize(250,20);
  _subFrame->AddFrame(caloHitLabel, lh1);
  _subFrame->AddFrame(caloHitBox, lh1);

  TGLabel *crvHitLabel  = new TGLabel(_subFrame, "CRV Hits");
  TGComboBox *crvHitBox = new TGComboBox(_subFrame,13);
  crvHitBox->Associate(this);
  crvHitBox->Resize(250,20);
  _subFrame->AddFrame(crvHitLabel, lh1);
  _subFrame->AddFrame(crvHitBox, lh1);

  TGLabel *trackLabel  = new TGLabel(_subFrame, "Tracks");
  TGListBox *trackBox = new TGListBox(_subFrame,12);
  trackBox->Associate(this);
  trackBox->Resize(250,60);
  trackBox->SetMultipleSelections(true);
  _subFrame->AddFrame(trackLabel, lh1);
  _subFrame->AddFrame(trackBox, lh1);

  _contentSelector=boost::shared_ptr<ContentSelector>(new ContentSelector(hitBox, caloHitBox, crvHitBox, trackBox,
                                                                          _g4ModuleLabel, _physicalVolumesMultiLabel,
                                                                          _protonBunchTimeLabel));

  TGHorizontalFrame *subFrameView1   = new TGHorizontalFrame(_subFrame,300,15);
  TGHorizontalFrame *subFrameView2   = new TGHorizontalFrame(_subFrame,300,15);
  TGTextButton *endViewButton       = new TGTextButton(subFrameView1, "End View", 70, buttoncontext);
  TGTextButton *sideViewButton      = new TGTextButton(subFrameView1, "Side View", 71, buttoncontext);
  TGTextButton *topViewButton       = new TGTextButton(subFrameView1, "Top View", 72, buttoncontext);
  TGTextButton *allTracksViewButton = new TGTextButton(subFrameView2, "All Tracks View", 73, buttoncontext);
  TGTextButton *resetViewButton     = new TGTextButton(subFrameView2, "Reset View", 74, buttoncontext);
  subFrameView1->AddFrame(endViewButton, lh1);
  subFrameView1->AddFrame(sideViewButton, lh1);
  subFrameView1->AddFrame(topViewButton, lh1);
  subFrameView2->AddFrame(allTracksViewButton, lh1);
  subFrameView2->AddFrame(resetViewButton, lh1);
  _subFrame->AddFrame(subFrameView1, lh0);
  _subFrame->AddFrame(subFrameView2, lh0);
  endViewButton->Associate(this);
  sideViewButton->Associate(this);
  topViewButton->Associate(this);
  allTracksViewButton->Associate(this);
  resetViewButton->Associate(this);

  TGHorizontalFrame *subFrameAnim = new TGHorizontalFrame(_subFrame,300,15);
  TGTextButton *animButtonStart   = new TGTextButton(subFrameAnim, "Start Animation", 40);
  TGTextButton *animButtonStop    = new TGTextButton(subFrameAnim, "Stop Animation",41);
  TGTextButton *animButtonReset   = new TGTextButton(subFrameAnim, "Reset",42);
  subFrameAnim->AddFrame(animButtonStart, lh1);
  subFrameAnim->AddFrame(animButtonStop, lh1);
  subFrameAnim->AddFrame(animButtonReset, lh1);
  _subFrame->AddFrame(subFrameAnim, lh0);
  animButtonStart->Associate(this);
  animButtonStop->Associate(this);
  animButtonReset->Associate(this);

  TGHorizontalFrame *subFrameAnimTime = new TGHorizontalFrame(_subFrame,300,15);
  TGLabel *timeIntervalLabel1  = new TGLabel(subFrameAnimTime, "Time Interval from");
  TGLabel *timeIntervalLabel2  = new TGLabel(subFrameAnimTime, "ns to");
  TGLabel *timeIntervalLabel3  = new TGLabel(subFrameAnimTime, "ns");
  _timeIntervalField1 = new TGTextEntry(subFrameAnimTime, new TGTextBuffer, 45);
  _timeIntervalField2 = new TGTextEntry(subFrameAnimTime, new TGTextBuffer, 46);
  subFrameAnimTime->AddFrame(timeIntervalLabel1, lh1);
  subFrameAnimTime->AddFrame(_timeIntervalField1, lh1);
  subFrameAnimTime->AddFrame(timeIntervalLabel2, lh1);
  subFrameAnimTime->AddFrame(_timeIntervalField2, lh1);
  subFrameAnimTime->AddFrame(timeIntervalLabel3, lh1);
  _subFrame->AddFrame(subFrameAnimTime, lh0);
  _timeIntervalField1->Associate(this);
  _timeIntervalField2->Associate(this);
  _timeIntervalField1->SetWidth(50);
  _timeIntervalField2->SetWidth(50);

  TGVerticalFrame *subFrameTimeWindow = new TGVerticalFrame(_subFrame,300,15);
  TGTextButton *allHitsTimeButton     = new TGTextButton(subFrameTimeWindow, "Time Window for all Hits", 80, buttoncontext);
  TGTextButton *allTracksTimeButton   = new TGTextButton(subFrameTimeWindow, "Time Window for all Tracks", 81, buttoncontext);
  subFrameTimeWindow->AddFrame(allHitsTimeButton, lh1);
  subFrameTimeWindow->AddFrame(allTracksTimeButton, lh1);
  _subFrame->AddFrame(subFrameTimeWindow, lh0);
  allHitsTimeButton->Associate(this);
  allTracksTimeButton->Associate(this);

  _repeatAnimationButton = new TGCheckButton(_subFrame,"Repeat Animation",43);
  _subFrame->AddFrame(_repeatAnimationButton, lh1);
  _repeatAnimationButton->Associate(this);

  TGHorizontalFrame *subFrameSave = new TGHorizontalFrame(_subFrame,300,15);
  TGTextButton *saveButton        = new TGTextButton(subFrameSave, "Save", 50);
  TGTextButton *saveAnimButton    = new TGTextButton(subFrameSave, "Save Animation", 51);
  subFrameSave->AddFrame(saveButton, lh1);
  subFrameSave->AddFrame(saveAnimButton, lh1);
  _subFrame->AddFrame(subFrameSave, lh0);
  saveButton->Associate(this);
  saveAnimButton->Associate(this);

  TGHorizontalFrame *subFrameSaveRootFile = new TGHorizontalFrame(_subFrame,300,15);
  TGTextButton *setRootFileButton         = new TGTextButton(subFrameSaveRootFile, "Set Root Tree File", 52);
  TGTextButton *addToRootFileButton       = new TGTextButton(subFrameSaveRootFile, "Add Event to Root Tree", 53);
  subFrameSaveRootFile->AddFrame(setRootFileButton, lh1);
  subFrameSaveRootFile->AddFrame(addToRootFileButton, lh1);
  _subFrame->AddFrame(subFrameSaveRootFile, lh0);
  setRootFileButton->Associate(this);
  addToRootFileButton->Associate(this);

  TGHorizontalFrame *subFrameFilterSetup = new TGHorizontalFrame(_subFrame,300,15);
  TGTextButton *filterButton = new TGTextButton(subFrameFilterSetup, "Filter", 63);
  TGTextButton *setupButton  = new TGTextButton(subFrameFilterSetup, "Setup", 64);
  subFrameFilterSetup->AddFrame(filterButton, lh1);
  subFrameFilterSetup->AddFrame(setupButton, lh1);
  _subFrame->AddFrame(subFrameFilterSetup, lh0);
  filterButton->Associate(this);
  setupButton->Associate(this);

  _eventInfo = new TGLabel*[2];
  for(int i=0; i<2; i++)
  {
    _eventInfo[i] = new TGLabel(_subFrame, "                      ");
    _eventInfo[i]->SetTextJustify(kTextLeft);
    _subFrame->AddFrame(_eventInfo[i], new TGLayoutHints(kLHintsLeft,2,0,2,1));
  }

  _footLine   = new TGHorizontalFrame(this,GetWidth(),GetHeight()*0.3);

//The following lines are a variation of something I found in CaloCellTimeMonitoring.C by Christophe Le Maner.
//Its purpose is to create a TRootEmbeddedCanvas with scrollbars.
  _infoCanvas = new TGCanvas(_footLine, GetWidth()-600, GetHeight()*0.28);
  _footLine->AddFrame(_infoCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,2,2,0,2));
  TGCompositeFrame *container = new TGCompositeFrame(_infoCanvas->GetViewPort());
  _infoCanvas->SetContainer(container);
  _infoEmbeddedCanvas = new TRootEmbeddedCanvas(0, container, 100, 100);//default size - exact size will be set later
  container->AddFrame(_infoEmbeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

  TGGroupFrame *zoomangleFrame  = new TGGroupFrame(_footLine,"Zoom & Angle");
  TGHorizontalFrame *zoomFrame1 = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGHorizontalFrame *zoomFrame2 = new TGHorizontalFrame(zoomangleFrame,500,50);
  TGHorizontalFrame *zoomFrame3 = new TGHorizontalFrame(zoomangleFrame,500,50);
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
  TGLabel *angleLabel2 = new TGLabel(angleFrame, "deg  theta");
  TGLabel *angleLabel3 = new TGLabel(angleFrame, "deg  psi");
  TGLabel *angleLabel4 = new TGLabel(angleFrame, "deg");
  _minXField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1501);
  _minYField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1502);
  _minZField = new TGTextEntry(zoomFrame1, new TGTextBuffer, 1503);
  _maxXField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1504);
  _maxYField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1505);
  _maxZField = new TGTextEntry(zoomFrame2, new TGTextBuffer, 1506);
  //_phiField   = new TGNumberEntry(angleFrame, 180., 3, 1601, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 360);
  //_thetaField   = new TGNumberEntry(angleFrame, 60., 3, 1602, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 360);
  //_psiField   = new TGNumberEntry(angleFrame, 90., 3, 1603, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 360);
  _phiField   = new TGNumberEntry(angleFrame, 180., 3, 1601, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -360, 640);
  _thetaField   = new TGNumberEntry(angleFrame, 60., 3, 1602, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -360, 640);
  _psiField   = new TGNumberEntry(angleFrame, 90., 3, 1603, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, -360, 640);
  _perspectiveButton = new TGRadioButton(perspectiveFrame, "perspective", 1700);
  _parallelButton    = new TGRadioButton(perspectiveFrame, "parallel", 1701);
  _zoomInButton      = new TGTextButton(zoomFrame3, "ZoomIn", 1507);
  _zoomOutButton     = new TGTextButton(zoomFrame3, "ZoomOut", 1508);
  _plusXButton       = new TGTextButton(zoomFrame3, "+X", 1509);
  _minusXButton      = new TGTextButton(zoomFrame3, "-X", 1510);
  _plusYButton       = new TGTextButton(zoomFrame3, "+Y", 1511);
  _minusYButton      = new TGTextButton(zoomFrame3, "-Y", 1512);
  _plusZButton       = new TGTextButton(zoomFrame3, "+Z", 1513);
  _minusZButton      = new TGTextButton(zoomFrame3, "-Z", 1514);
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
  zoomFrame3->AddFrame(_zoomInButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_zoomOutButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_plusXButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_minusXButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_plusYButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_minusYButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_plusZButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));
  zoomFrame3->AddFrame(_minusZButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,1,0,1,0));

  zoomangleFrame->AddFrame(zoomFrame1, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(zoomFrame2, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(setRangeButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(angleFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(setAngleButton, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(perspectiveFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  zoomangleFrame->AddFrame(zoomFrame3, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  _footLine->AddFrame(zoomangleFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));

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
  _zoomInButton->Associate(this);
  _zoomOutButton->Associate(this);
  _plusXButton->Associate(this);
  _minusXButton->Associate(this);
  _plusYButton->Associate(this);
  _minusYButton->Associate(this);
  _plusZButton->Associate(this);
  _minusZButton->Associate(this);
  _perspectiveButton->SetState(kButtonDown);
  _parallelButton->SetState(kButtonUp);

  TGVerticalFrame *innerFrame1   = new TGVerticalFrame(_footLine,100,400);

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
  TGTextButton *updateButton       = new TGTextButton(navigationFrame, "&Update", 1112);
  navigationFrame->AddFrame(quitButton, new TGLayoutHints(kLHintsLeft,3,0,3,0));
  navigationFrame->AddFrame(nextButton, new TGLayoutHints(kLHintsLeft,3,0,3,0));
  navigationFrame->AddFrame(updateButton, new TGLayoutHints(kLHintsLeft,3,0,3,0));
  innerFrame1->AddFrame(navigationFrame, new TGLayoutHints(kLHintsLeft,3,0,3,0));

  quitButton->Associate(this);
  nextButton->Associate(this);
  updateButton->Associate(this);
  _minHitField->Associate(this);
  _eventToFindField->Associate(this);
  applyButton->Associate(this);
  goButton->Associate(this);

  mu2e::ConfigFileLookupPolicy configFile;
  std::string logoFileName = configFile("Offline/EventDisplay/src/logo_small.png");
  const TGPicture *logo = gClient->GetPicture(logoFileName.c_str());
  TGIcon *icon = new TGIcon(navigationFrame, logo, 50, 50);
  navigationFrame->AddFrame(icon, new TGLayoutHints(kLHintsLeft,20,0,0,0));

  _footLine->AddFrame(innerFrame1, new TGLayoutHints(kLHintsLeft,0,0,0,0));
  AddFrame(_footLine, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  MapSubwindows();
  SetWindowName("Mu2e Event Display");
  MapWindow();

  _mainCanvas->GetCanvas()->cd();
  _mainPad = new TPad("mainPad","Detector", 0, 0, 1, 1, 5,1,1);
  _mainPad->SetFillColor(1);
  _mainPad->Draw();

  _infoEmbeddedCanvas->GetCanvas()->cd();
  _infoPad = new TPad("infoPad","InfoField", 0, 0, 1, 1, 5,1,1);
  _infoPad->SetFillColor(0);
  _infoPad->Draw();

  for(int i=0; i<20; i++)
  {
    float r,g,b;
    TColor::HLS2RGB(i*360/20,.5,.5,r,g,b);
    if(!gROOT->GetColor(i+2000)) new TColor(i+2000,r,g,b);
  }

  _mainPad->cd();
  _dataInterface = boost::shared_ptr<DataInterface>(new DataInterface(this,_kalStepSize));
  _rootFileManager = boost::shared_ptr<RootFileManager>(new RootFileManager);
  _rootFileManagerAnim = boost::shared_ptr<RootFileManager>(new RootFileManager);

//syntax for Format from http://root.cern.ch/phpBB3/viewtopic.php?t=8700
  gPad->AddExec("keyboardInput",TString::Format("((mu2e_eventdisplay::EventDisplayFrame*)%p)->keyboardInput()",static_cast<void*>(this)));
}

void EventDisplayFrame::initSetup()
{
  _whiteBackground=false;
  _useHitColors=true;
  _useTrackColors=true;
  _showSupportStructures=true;
  _showCRV=false;
  _showOtherStructures=false;
  _showMuonBeamStop=false;
  _showProtonAbsorber=false;
  if(_wideband || _extracted)
  {
    _showSupportStructures=false;
    _showCRV=true;
  }
}

void EventDisplayFrame::changeSetup(bool whiteBackground, bool useHitColors, bool useTrackColors,
                                    bool showSupportStructures,
                                    bool showCRV,
                                    bool showOtherStructures,
                                    bool showMuonBeamStop,
                                    bool showProtonAbsorber)
{
  _mainPad->cd();
  if(_whiteBackground!=whiteBackground)
  {
    _whiteBackground=whiteBackground;
    if(whiteBackground) _mainPad->SetFillColor(0);
    else _mainPad->SetFillColor(1);
    _dataInterface->useHitColors(useHitColors, whiteBackground);
    _dataInterface->useTrackColors(_contentSelector, useTrackColors, whiteBackground);
    updateHitLegend(useHitColors);
    updateTrackLegend(useTrackColors);
  }

  if(_useHitColors!=useHitColors)
  {
    _useHitColors=useHitColors;
    _dataInterface->useHitColors(useHitColors, whiteBackground);
    updateHitLegend(useHitColors);
  }

  if(_useTrackColors!=useTrackColors)
  {
    _useTrackColors=useTrackColors;
    _dataInterface->useTrackColors(_contentSelector, useTrackColors, whiteBackground);
    updateTrackLegend(useTrackColors);
  }

  bool redraw=false;
  if(_showSupportStructures!=showSupportStructures)
  {
     redraw=true;
     _showSupportStructures=showSupportStructures;
     _dataInterface->makeSupportStructuresVisible(showSupportStructures);
  }
  if(_showCRV!=showCRV)
  {
     redraw=true;
     _showCRV=showCRV;
     _dataInterface->makeCrvScintillatorBarsVisible(showCRV);
  }
  if(_showOtherStructures!=showOtherStructures)
  {
     redraw=true;
     _showOtherStructures=showOtherStructures;
     _dataInterface->makeOtherStructuresVisible(showOtherStructures);
  }
  if(_showMuonBeamStop!=showMuonBeamStop)
  {
     redraw=true;
     _showMuonBeamStop=showMuonBeamStop;
     _dataInterface->makeMuonBeamStopStructuresVisible(showMuonBeamStop);
  }
  if(_showProtonAbsorber!=showProtonAbsorber)
  {
     redraw=true;
     _showProtonAbsorber=showProtonAbsorber;
     _dataInterface->makeProtonAbsorberVisible(showProtonAbsorber);
  }

  if(std::isnan(_timeCurrent) || redraw) drawEverything();
  else drawSituation();
}

EventDisplayFrame::~EventDisplayFrame()
{
  // TODO
  // delete timer;
  // Cleanup();
}

void EventDisplayFrame::keyboardInput()
{
  EventDisplayViewSetup::input();
  fillZoomAngleFields();
}

Bool_t EventDisplayFrame::HandleConfigureNotify(Event_t *event)
{
// This is a modified version of the function from TGFrame.cxx
   if ((event->fWidth != fWidth) || (event->fHeight != fHeight))
   {
      fWidth  = event->fWidth;
      fHeight = event->fHeight;

      _mainFrame->SetWidth(fWidth-270);
      _mainFrame->SetHeight(fHeight*0.7);
      _subFrame->SetWidth(270);
      _subFrame->SetHeight(fHeight*0.7);
      _footLine->SetWidth(fWidth);
      _footLine->SetHeight(fHeight*0.3);

      _mainCanvas->SetWidth(fWidth-270);
      _mainCanvas->SetHeight(fHeight*0.7);
      _infoCanvas->SetWidth(fWidth-600);
      _infoCanvas->SetHeight(fHeight*0.28);

      Layout();
   }
   return kTRUE;
}

void EventDisplayFrame::fillZoomAngleFields()
{
  if(_mainPad->GetView()==nullptr) return;
  double min[3], max[3];
  _mainPad->GetView()->GetRange(min,max);
  _minXField->SetText(Form("%.0f",min[0]));
  _minYField->SetText(Form("%.0f",min[1]));
  _minZField->SetText(Form("%.0f",min[2]));
  _maxXField->SetText(Form("%.0f",max[0]));
  _maxYField->SetText(Form("%.0f",max[1]));
  _maxZField->SetText(Form("%.0f",max[2]));
  _phiField->SetText(Form("%.0f",_mainPad->GetView()->GetLongitude()));
  _thetaField->SetText(Form("%.0f",_mainPad->GetView()->GetLatitude()));
  _psiField->SetText(Form("%.0f",_mainPad->GetView()->GetPsi()));
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
  DataInterface::spaceminmax m=_dataInterface->getSpaceBoundary(true, true, true, false, _showCRV);
  if(_wideband) m=_dataInterface->getSpaceBoundary(false, false, false, false, true);
  if(_extracted) m=_dataInterface->getSpaceBoundary(true, false, true, false, true);
  _mainPad->GetView()->SetRange(m.minx,m.miny,m.minz,m.maxx,m.maxy,m.maxz);
  if(_wideband || _extracted)
  {
    EventDisplayViewSetup::sideview();
    _parallelButton->SetState(kButtonDown);
    _perspectiveButton->SetState(kButtonUp);
    _mainPad->GetView()->SetParallel();
  }
  else EventDisplayViewSetup::perspectiveview();
  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::setEvent(const art::Event& event, bool firstLoop, const mu2e::CRVCalib &calib)
{

  _eventNumber=event.id().event();
  _subrunNumber=event.id().subRun();
  _runNumber=event.id().run();

  _dataInterface->setCRVCalib(calib);

  _contentSelector->setAvailableCollections(event);
  if(firstLoop) _contentSelector->firstLoop();
  fillEvent(firstLoop);
  fillZoomAngleFields();

  gApplication->Run(true);
}

void EventDisplayFrame::fillEvent(bool firstLoop)
{
  _findEvent=false;
  _mainPad->cd();

  std::string eventInfoText;
  eventInfoText=Form("Event #: %i",_eventNumber);
  if(_eventNumberText==nullptr)
  {
    _eventNumberText = new TText(0.6,-0.8,eventInfoText.c_str());
    _eventNumberText->SetTextColor(5);
    _eventNumberText->SetTextSize(0.025);
    _eventNumberText->Draw("same");
  }
  else _eventNumberText->SetTitle(eventInfoText.c_str());
  eventInfoText=Form("Sub Run #: %i",_subrunNumber);
  if(_subrunNumberText==nullptr)
  {
    _subrunNumberText = new TText(0.6,-0.75,eventInfoText.c_str());
    _subrunNumberText->SetTextColor(5);
    _subrunNumberText->SetTextSize(0.025);
    _subrunNumberText->Draw("same");
  }
  else _subrunNumberText->SetTitle(eventInfoText.c_str());
  eventInfoText=Form("Run #: %i",_runNumber);
  if(_runNumberText==nullptr)
  {
    _runNumberText = new TText(0.6,-0.7,eventInfoText.c_str());
    _runNumberText->SetTextColor(5);
    _runNumberText->SetTextSize(0.025);
    _runNumberText->Draw("same");
  }
  else _runNumberText->SetTitle(eventInfoText.c_str());

  _dataInterface->fillEvent(_contentSelector);
  _dataInterface->useHitColors(_useHitColors, _whiteBackground);
  _dataInterface->useTrackColors(_contentSelector, _useTrackColors, _whiteBackground);
  _dataInterface->makeCrvScintillatorBarsVisible(_showCRV);

  updateTimeIntervalFields();
  updateHitLegend(_useHitColors);
  updateTrackLegend(_useTrackColors);

  _eventInfo[0]->SetText(Form("Number of tracker hits: %i",_dataInterface->getNumberHits()));
  _eventInfo[1]->SetText(Form("Number of calorimeter hits: %i",_dataInterface->getNumberCrystalHits()));
  this->Layout();

  drawEverything();
}

void EventDisplayFrame::updateTimeIntervalFields(bool allTracks)
{
  double mint, maxt;
  if(allTracks)
  {
    mint=_dataInterface->getTracksTimeBoundary().mint;
    maxt=_dataInterface->getTracksTimeBoundary().maxt;
  }
  else
  {
    mint=_dataInterface->getHitsTimeBoundary().mint;
    maxt=_dataInterface->getHitsTimeBoundary().maxt;
  }

  if(maxt<mint) {mint=NAN; maxt=NAN;}
  _timeIntervalField1->SetText(Form("%.0f",mint));
  _timeIntervalField2->SetText(Form("%.0f",maxt));
}

void EventDisplayFrame::updateHitLegend(bool draw)
{
  for(int i=0; i<21; i++)
  {
    delete _legendBox[i];
    delete _legendText[i];
    _legendBox[i]=nullptr;
    _legendText[i]=nullptr;
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
        if(maxt<=mint) t=NAN;
        _legendText[i]=new TText(0.72,0.54+i*0.02,Form("%+.3e ns",t));
        _legendText[i]->SetTextColor(_whiteBackground?kBlack:kWhite);
        _legendText[i]->SetTextSize(0.025);
        _legendText[i]->Draw();
      }
    }
  }
}

void EventDisplayFrame::updateTrackLegend(bool draw)
{
  for(int i=0; i<30; i++)
  {
    delete _legendParticleGroup[i];
    delete _legendParticleLine[i];
    delete _legendParticleText[i];
    _legendParticleGroup[i]=nullptr;
    _legendParticleLine[i]=nullptr;
    _legendParticleText[i]=nullptr;
  }

  if(draw)
  {
    std::vector<ContentSelector::trackInfoStruct> selectedTracks=_contentSelector->getSelectedTrackNames();
    TrackColorSelector colorSelector(&selectedTracks, _whiteBackground);
    colorSelector.drawTrackLegend(_legendParticleGroup, _legendParticleText, _legendParticleLine);
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
    case kC_TEXTENTRY:
      switch (GET_SUBMSG(msg))
      {
        case kTE_TEXTCHANGED:
          if (param1 ==1601 || param1 == 1602 || param1 == 1603)
          {
            double phi = _phiField->GetNumber();
            double theta = _thetaField->GetNumber();
            double psi = _psiField->GetNumber();
            while (phi>360)
            {
              phi-=360;
              _phiField->SetText(Form("%.0f",phi));
            }
            while (phi<0)
            {
              phi+=360;
              _phiField->SetText(Form("%.0f",phi));
            }
            while (theta>360)
            {
              theta-=360;
              _thetaField->SetText(Form("%.0f",theta));
            }
            while (theta<0)
            {
              theta+=360;
              _thetaField->SetText(Form("%.0f",theta));
            }
            while (psi>360) {
              psi-=360;
              _psiField->SetText(Form("%.0f",psi));
            }
            while (psi<0)
            {
              psi+=360;
              _psiField->SetText(Form("%.0f",psi));
            }

            int irep=0;
            _mainPad->GetView()->SetView(phi,theta,psi,irep);
            _mainPad->SetPhi(-90-phi);
            _mainPad->SetTheta(90-theta);
            _mainPad->Modified();
            _mainPad->Update();
          }
          break;
        default:
          break;
      }
      break;

    case kC_COMMAND:
      switch (GET_SUBMSG(msg))
      {
        case kCM_BUTTON: if(param1==1001) CloseWindow();
                         if(param1==40) prepareAnimation();
                         if(param1==41)
                         {
                           _timer->Stop();
                           _timeCurrent=NAN;
                           if(_saveAnim) createAnimFile();
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
                           gApplication->Terminate();
                         }
                         if(param1==1112)
                         {
                           _timer->Stop();
                           _timeCurrent=NAN;
                           fillEvent();
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
                         if(param1>=70 && param1<=74)
                         {
                           _mainPad->cd();
                           DataInterface::spaceminmax m=_dataInterface->getSpaceBoundary(true, true, true, param1==73, _showCRV);
                           if(_wideband) m=_dataInterface->getSpaceBoundary(false, false, false, false, true);
                           if(_extracted) m=_dataInterface->getSpaceBoundary(true, false, true, false, true);
                           _mainPad->GetView()->SetRange(m.minx,m.miny,m.minz,m.maxx,m.maxy,m.maxz);
                           if(param1<73)
                           {
                             if(param1==70) EventDisplayViewSetup::endview();
                             if(param1==71) EventDisplayViewSetup::sideview();
                             if(param1==72) EventDisplayViewSetup::topview();
                             _parallelButton->SetState(kButtonDown);
                             _perspectiveButton->SetState(kButtonUp);
                             _mainPad->GetView()->SetParallel();
                           }
                           if(param1==74)
                           {
                             EventDisplayViewSetup::perspectiveview();
                             _parallelButton->SetState(kButtonUp);
                             _perspectiveButton->SetState(kButtonDown);
                             _mainPad->GetView()->SetPerspective();
                           }
                           _mainPad->Modified();
                           _mainPad->Update();
                           fillZoomAngleFields();
                         }
                         if(param1==80) updateTimeIntervalFields(false);
                         if(param1==81) updateTimeIntervalFields(true);
                         if(param1==50)
                         {
                           std::string filename;
                           if(SaveDialogManager::singleImage(filename)) _mainPad->GetCanvas()->SaveAs(filename.c_str());
                         }
                         if(param1==51)
                         {
                           std::string filename;
                           bool isRootFile=false;
                           if(SaveDialogManager::animatedImage(filename,isRootFile))
                           {
                             _saveAnim=true;
                             if(isRootFile) _saveAnimRoot=true; else _saveAnimRoot=false;
                             _saveAnimFile=filename;
                             prepareAnimation();
                           }
                         }
                         if(param1==52)
                         {
                           std::string filename;
                           if(SaveDialogManager::rootTree(filename)) _rootFileManager->setFile(filename.c_str());
                         }
                         if(param1==53) _rootFileManager->storeEvent(_mainPad);
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
                           double phi = _phiField->GetNumber();
                           double theta = _thetaField->GetNumber();
                           double psi = _psiField->GetNumber();
                           int irep=0;
                           _mainPad->GetView()->SetView(phi,theta,psi,irep);
                           _mainPad->SetPhi(-90.0-phi);
                           _mainPad->SetTheta(90.0-theta);
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
                         if(param1==1507 || param1==1508)
                         {
                           double min[3],max[3],hdist,cen;
                           _mainPad->GetView()->GetRange(min,max);
                           if(min[0]<max[0] && min[1]<max[1] && min[2]<max[2])
                           {
                             for(int i=0; i<3; i++)
                             {
                               hdist = (max[i] - min[i])/2.0;
                               cen   = min[i] + hdist;
                               if(param1==1507) hdist*=0.9;
                               else hdist/=0.9;
                               min[i] = cen - hdist;
                               max[i] = cen + hdist;
                             }
                             _mainPad->GetView()->SetRange(min,max);
                             _minXField->SetText(Form("%.0f",min[0]));
                             _minYField->SetText(Form("%.0f",min[1]));
                             _minZField->SetText(Form("%.0f",min[2]));
                             _maxXField->SetText(Form("%.0f",max[0]));
                             _maxYField->SetText(Form("%.0f",max[1]));
                             _maxZField->SetText(Form("%.0f",max[2]));
                             _mainPad->Modified();
                             _mainPad->Update();
                           }
                         }
                         if(1509 <= param1 && param1 <= 1514)
                         {
                           double min[3],max[3];
                           min[0]=atof(_minXField->GetText());
                           min[1]=atof(_minYField->GetText());
                           min[2]=atof(_minZField->GetText());
                           max[0]=atof(_maxXField->GetText());
                           max[1]=atof(_maxYField->GetText());
                           max[2]=atof(_maxZField->GetText());
                           int mode = (param1 - 1509)/2;
                           int dir = (param1 - 1509)%2;
                           double diff = 100;
                           if (!dir) diff*=-1;
                           if (0 <= mode && mode <3)
                           {
                             min[mode] += diff;
                             max[mode] += diff;
                             _mainPad->GetView()->SetRange(min,max);
                             if(mode == 0)
                             {
                               _minXField->SetText(Form("%.0f",min[0]));
                               _maxXField->SetText(Form("%.0f",max[0]));
                             }
                             else if(mode == 1)
                             {
                               _minYField->SetText(Form("%.0f",min[1]));
                               _maxYField->SetText(Form("%.0f",max[1]));
                             }
                             else if(mode == 2)
                             {
                               _minZField->SetText(Form("%.0f",min[2]));
                               _maxZField->SetText(Form("%.0f",max[2]));
                             }
                             else {}
                             _mainPad->Modified();
                             _mainPad->Update();
                           }
                         }
                         if(param1==63)
                         {
                           new FilterDialog(gClient->GetRoot(), _dataInterface, _contentSelector);
                           _dataInterface->useHitColors(_useHitColors, _whiteBackground);
                           updateHitLegend(_useHitColors);
                           drawEverything();
                         }
                         if(param1==64)
                         {
                           new SetupDialog(gClient->GetRoot(), this, _whiteBackground,
                                           _useHitColors, _useTrackColors,
                                           _showSupportStructures, _showCRV, _showOtherStructures,
                                           _showMuonBeamStop, _showProtonAbsorber);
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
/*
  case kCM_COMBOBOX : if(param1==10) fillEvent();
                      if(param1==11) fillEvent();
                      if(param1==13) fillEvent();
                      break;
  case kCM_LISTBOX : if(param1==12) fillEvent();
                     break;
*/
      }
      break;
  }
  return kTRUE;
}

void EventDisplayFrame::prepareAnimation()
{
  _timer->Stop();           //needed if an animation is already running
  delete _clock; _clock=nullptr;
  _mainPad->cd();
  _dataInterface->startComponents();
  _mainPad->Modified();
  _mainPad->Update();
  _timeStart=atof(_timeIntervalField1->GetText());
  _timeStop=atof(_timeIntervalField2->GetText());

  if(std::isnan(_timeStart) || std::isnan(_timeStop) || (_timeStop-_timeStart)<=0.0) return;

  if(_saveAnim && _saveAnimRoot) _rootFileManagerAnim->setFile(_saveAnimFile.c_str());
  double diff=_timeStop-_timeStart;
  _timeStart-=diff*0.01;
  _timeStop+=diff*0.01;
  _timeCurrent=_timeStart;
  _saveAnimCounter=0;
  _timer->Start(100,kFALSE);
}


Bool_t EventDisplayFrame::HandleTimer(TTimer *)
{
  _timer->Reset();
  _timeCurrent+=(_timeStop-_timeStart)/50;
  drawSituation();
  if(_saveAnim)
  {
    _saveAnimCounter++;
    if(_saveAnimRoot) _rootFileManagerAnim->storeEvent(_mainPad);
    else
    {
//      if(_saveAnimCounter%3==0) //save only every 3rd gif to make final file smaller
      {
        _mainPad->GetCanvas()->SaveAs(Form("%s.tmp_%04i.gif",_saveAnimFile.c_str(),_saveAnimCounter));
      }
    }
  }
  if(_timeCurrent>=_timeStop)
  {
    _timer->Stop();
    _timeCurrent=NAN;
    if(_saveAnim) createAnimFile();
    if(_repeatAnimationButton->GetState()==kButtonDown) prepareAnimation();
  }
  return kTRUE;
}

void EventDisplayFrame::createAnimFile()
{
  _saveAnim=false;
  if(_saveAnimRoot) _rootFileManagerAnim->write();
  else
  {
    gSystem->Exec(Form("convert -delay 50 -loop 0 %s.tmp_*.gif %s",_saveAnimFile.c_str(),_saveAnimFile.c_str()));
    gSystem->Exec(Form("rm -f %s.tmp_*.gif",_saveAnimFile.c_str()));
  }
}

void EventDisplayFrame::drawSituation()
{
  _mainPad->cd();
  _dataInterface->updateComponents(_timeCurrent, _contentSelector);
  if(!_clock)
  {
    _clock = new TText(0.6,-0.9,Form("%+.4e ns",_timeCurrent));
    _clock->SetTextColor(5);
    _clock->SetTextSize(0.025);
    _clock->Draw("same");
  }
  else
  {
    _clock->SetTitle(Form("%+.4e ns",_timeCurrent));
  }

  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::drawEverything()
{
  _mainPad->cd();
  delete _clock; _clock=nullptr;
  if(TAxis3D::GetPadAxis(_mainPad)==nullptr) _mainPad->GetView()->ShowAxis();
  TAxis3D::GetPadAxis(_mainPad)->SetLabelSize(0.025);
  _mainPad->Modified();
  _mainPad->Update();
  _dataInterface->updateComponents(NAN, _contentSelector);
  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::showInfo(TObject *o)  //ROOT accepts only bare pointers here
{
  ComponentInfoContainer *container=dynamic_cast<ComponentInfoContainer*>(o);
  if(!container) return;

  _infoPad->cd();
  _infoPad->Clear();
  unsigned int w,h;
  container->getComponentInfo()->getExpectedSize(w,h);
  if(w<_infoCanvas->GetWidth()-20) w=_infoCanvas->GetWidth()-20;
  if(h<_infoCanvas->GetHeight()-20) h=_infoCanvas->GetHeight()-20;
  if(w>10000) w=10000;
  if(h>10000) h=10000;
  _infoEmbeddedCanvas->SetWidth(w);
  _infoEmbeddedCanvas->SetHeight(h);
  _infoCanvas->Layout();
  _infoPad->cd();
  _infoPad->Range(0, 0, 1, 1);
  _infoPad->Modified();
  _infoPad->Update();
  container->getComponentInfo()->showInfo(w,h);
  _infoPad->Modified();
  _infoPad->Update();
  _mainPad->cd();
}

void EventDisplayFrame::addHistDraw()
{
  unsigned int i = _histDrawVector.size();
  boost::shared_ptr<HistDraw> histDraw(new HistDraw(i,_infoPad,_mainPad,_infoEmbeddedCanvas,_infoCanvas));
  _histDrawVector.push_back(histDraw);
}

}
