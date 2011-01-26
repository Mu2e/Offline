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

#include "EventDisplayFrame.h"
#include "VirtualShape.h"
#include "DataInterface.h"

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
  MoveResize(20,20,width-50,height-50);

  _timer=new TTimer();
  _timer->SetObject(this);
  _clock=NULL;
  _isClosed=false;
  _saveAnim=false;
  for(int i=0; i<30; i++)
  {
    _legendText[i]=NULL;
    _legendBox[i]=NULL;
  }

  //bare pointers needed since ROOT manages the following object
  TGHorizontalFrame *mainFrame = new TGHorizontalFrame(this,800,400);
  _mainCanvas = new TRootEmbeddedCanvas("EventDisplayCanvas",mainFrame,GetWidth()-300,GetHeight()-150);
  TGVerticalFrame *subFrame = new TGVerticalFrame(mainFrame,300,400);

  mainFrame->AddFrame(_mainCanvas, new TGLayoutHints(kLHintsTop));
  mainFrame->AddFrame(subFrame, new TGLayoutHints(kLHintsTop));
  AddFrame(mainFrame, new TGLayoutHints(kLHintsTop, 2,2,2,2));

  TGLayoutHints *lh0 = new TGLayoutHints(kLHintsTop,0,0,0,0);
  TGLayoutHints *lh1 = new TGLayoutHints(kLHintsTop,2,1,2,2);

  TGHorizontalFrame *hitFrame  = new TGHorizontalFrame(subFrame,200,15);
  TGLabel *hitLabel  = new TGLabel(hitFrame, "hits");
  TGComboBox *hitBox = new TGComboBox(hitFrame,10);
  hitBox->Associate(this);
  hitBox->AddEntry("",0);
  hitBox->Select(0);
  hitBox->Resize(150,20);
  hitFrame->AddFrame(hitLabel, lh1);
  hitFrame->AddFrame(hitBox, lh1);
  subFrame->AddFrame(hitFrame, lh0);

  _unhitButton = new TGCheckButton(subFrame,"Show Unhit Straws",31);
  subFrame->AddFrame(_unhitButton, lh1);
  _unhitButton->Associate(this);

  _supportStructureButton = new TGCheckButton(subFrame,"Show Support Structures",32);
  _supportStructureButton->SetState(kButtonDown);
  subFrame->AddFrame(_supportStructureButton, lh1);
  _supportStructureButton->Associate(this);

  _outsideTracksButton = new TGCheckButton(subFrame,"Show Tracks Outside of Detector",33);
  subFrame->AddFrame(_outsideTracksButton, lh1);
  _outsideTracksButton->Associate(this);

  TGListBox **particleBox   = new TGListBox*[2];
  for(int i=0; i<2; i++)
  {
    particleBox[i] = new TGListBox(subFrame,20+i);
    particleBox[i]->Associate(this);
    particleBox[i]->Resize(150,60);
    particleBox[i]->SetMultipleSelections(true);
    subFrame->AddFrame(particleBox[i], lh1);
  }

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

  TGCheckButton *loopButton = new TGCheckButton(subFrame,"Repeat Animation",42);
  subFrame->AddFrame(loopButton, lh1);
  loopButton->Associate(this);

  TGTextButton *view3dButton = new TGTextButton(subFrame, "View3D", 60);
  subFrame->AddFrame(view3dButton, lh1);
  view3dButton->Associate(this);

  TGHorizontalFrame *subFrameSave = new TGHorizontalFrame(subFrame,300,15);
  TGTextButton *saveButton        = new TGTextButton(subFrameSave, "Save", 50);
  TGTextButton *saveAnimButton    = new TGTextButton(subFrameSave, "Save Animation", 51);
  subFrameSave->AddFrame(saveButton, lh1);
  subFrameSave->AddFrame(saveAnimButton, lh1);
  subFrame->AddFrame(subFrameSave, lh0);
  saveButton->Associate(this);
  saveAnimButton->Associate(this);

  TGCheckButton *backgroundButton = new TGCheckButton(subFrame,"White Background",33);
  subFrame->AddFrame(backgroundButton, lh1);
  backgroundButton->Associate(this);
  backgroundButton->SetState(kButtonDown);
  _backgroundColor=0;

  _eventInfo = new TGLabel*[3];
  for(int i=0; i<3; i++)
  {
    _eventInfo[i] = new TGLabel(subFrame, "Place Holder for Event Info");
    _eventInfo[i]->SetTextJustify(kTextLeft);
    subFrame->AddFrame(_eventInfo[i], new TGLayoutHints(kLHintsLeft,2,0,2,1));
  }

  TGHorizontalFrame *footLine   = new TGHorizontalFrame(this,800,100);
  _infoCanvas = new TRootEmbeddedCanvas("InfoCanvas",footLine,GetWidth()-520,GetHeight()-430);
  footLine->AddFrame(_infoCanvas, new TGLayoutHints(kLHintsTop));

  TGVerticalFrame *navigationFrame = new TGVerticalFrame(footLine,100,200);
  TGTextButton *quitButton         = new TGTextButton(navigationFrame, "&Quit", 1001);
  TGTextButton *nextButton         = new TGTextButton(navigationFrame, "&Next", 1111);
  navigationFrame->AddFrame(quitButton, new TGLayoutHints(kLHintsLeft,10,0,10,0));
  navigationFrame->AddFrame(nextButton, new TGLayoutHints(kLHintsLeft,10,0,10,0));
  footLine->AddFrame(navigationFrame, new TGLayoutHints(kLHintsLeft,10,0,10,0));

  TGGroupFrame *optionsFrame     = new TGGroupFrame(footLine,"Options");
  TGHorizontalFrame *filterFrame = new TGHorizontalFrame(optionsFrame,500,50);
  TGLabel *minHitLabel     = new TGLabel(filterFrame, "minimum hits");
  _minHitField = new TGTextEntry(filterFrame, new TGTextBuffer, 1101);
  _minHits=0;  
  _minHitField->SetWidth(100); 
  _minHitField->SetText("0");
  TGTextButton *applyButton = new TGTextButton(filterFrame, "&Apply", 1100);
  filterFrame->AddFrame(minHitLabel, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,15,0,5,0));
  filterFrame->AddFrame(_minHitField, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,5,0,10,0));
  filterFrame->AddFrame(applyButton, new TGLayoutHints(kLHintsLeft|kLHintsCenterY,10,0,5,0));
  optionsFrame->AddFrame(filterFrame, new TGLayoutHints(kLHintsLeft,5,0,0,0));
  footLine->AddFrame(optionsFrame, new TGLayoutHints(kLHintsLeft,0,0,0,0));

  quitButton->Associate(this);
  nextButton->Associate(this);
  _minHitField->Associate(this);
  applyButton->Associate(this);

  std::string logoFileName=getenv("MU2E_BASE_RELEASE");
  logoFileName.append("/EventDisplay/src/logo.png");
  TGPicture *logo = (TGPicture *)gClient->GetPicture(logoFileName.c_str());
  TGIcon *icon = new TGIcon(footLine, logo, 50, 50);
  footLine->AddFrame(icon, new TGLayoutHints(kLHintsLeft,10,10,10,10));

  AddFrame(footLine, new TGLayoutHints(kLHintsLeft,1,1,1,1));

  MapSubwindows();
  SetWindowName("Mu2e Event Display");
  MapWindow();

  _mainCanvas->GetCanvas()->cd();
  _mainPad = new TPad("mainPad","Detector", 0, 0, 1, 1, 5,1,1);  
  _mainPad->SetFillColor(_backgroundColor);
  _mainPad->Draw();

  _infoCanvas->GetCanvas()->cd();
  _infoPad = new TPad("infoPad","InfoField", 0, 0, 1, 1, 5,1,1);  
  _infoPad->SetFillColor(_backgroundColor);
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
//  mainPad->SetView(new TView3D());  //not needed - gets provided by TGeoManager
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
      _mainCanvas->SetWidth(fWidth-300);
      _mainCanvas->SetHeight(fHeight-150);
      _infoCanvas->SetWidth(fWidth-520);
      _infoCanvas->SetHeight(fHeight-430);
      Layout();
   }
   return kTRUE;
}

void EventDisplayFrame::fillGeometry()
{
  _mainPad->cd();
  _dataInterface->fillGeometry();
  double minx=_dataInterface->getHitsBoundary().minx;
  double miny=_dataInterface->getHitsBoundary().miny;
  double minz=_dataInterface->getHitsBoundary().minz;
  double maxx=_dataInterface->getHitsBoundary().maxx;
  double maxy=_dataInterface->getHitsBoundary().maxy;
  double maxz=_dataInterface->getHitsBoundary().maxz;
  _mainPad->GetView()->SetRange(minx,miny,minz,maxx,maxy,maxz);
  _mainPad->GetView()->AdjustScales();
  _mainPad->Modified();
  _mainPad->Update();
}

void EventDisplayFrame::fillEvent(const edm::Event& event)
{
  _mainPad->cd();
  _dataInterface->fillEvent(event);
  if(_outsideTracksButton->GetState()==kButtonDown)
  {
    double minx=_dataInterface->getTracksBoundary().minx;
    double miny=_dataInterface->getTracksBoundary().miny;
    double minz=_dataInterface->getTracksBoundary().minz;
    double maxx=_dataInterface->getTracksBoundary().maxx;
    double maxy=_dataInterface->getTracksBoundary().maxy;
    double maxz=_dataInterface->getTracksBoundary().maxz;
    _mainPad->GetView()->SetRange(minx,miny,minz,maxx,maxy,maxz);
    _mainPad->GetView()->AdjustScales();;
  }

  for(int i=0; i<21; i++)
  {
    if(i<20)
    {
      if(_legendBox[i]!=NULL) delete _legendBox[i];
      _legendBox[i]=new TBox(0.6,0.55+i*0.02,0.7,0.57+i*0.02);
      _legendBox[i]->SetFillColor(i+2000);
      _legendBox[i]->Draw();
    }
    if(i%4==0)
    {
      double mint=_dataInterface->getHitsBoundary().mint;
      double maxt=_dataInterface->getHitsBoundary().maxt;
      double t=i*(maxt-mint)/20.0+mint;
      char s[50];
      sprintf(s,"%+.3e ns",t);
      if(_legendText[i]!=NULL) delete _legendText[i];
      _legendText[i]=new TText(0.72,0.54+i*0.02,s);
      _legendText[i]->SetTextColor(1);
      _legendText[i]->SetTextSize(0.025);
      _legendText[i]->Draw();
    }
  }

  char eventInfoText[50];
  sprintf(eventInfoText,"Event #: %i",event.id().event());
  _eventInfo[0]->SetText(eventInfoText);
  sprintf(eventInfoText,"Run #: %i",event.id().run());
  _eventInfo[1]->SetText(eventInfoText);
  sprintf(eventInfoText,"Number of hit straws: %i",_dataInterface->getNumberHits());
  _eventInfo[2]->SetText(eventInfoText);
  this->Layout();

  drawEverything();
  gApplication->Run(true);
}

bool EventDisplayFrame::isClosed() const
{
  return _isClosed;
}

int EventDisplayFrame::getMinimumHits() const
{
  return _minHits;
}

void EventDisplayFrame::CloseWindow()
{
  _isClosed=true;
  _timer->Stop();
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
                           if(_saveAnim) combineAnimFiles();
                         }
                         if(param1==42) 
                         {
                           _timer->Stop(); 
                           drawEverything();
                         }
                         if(param1==1111)
                         {
                           _timer->Stop();
                           gApplication->Terminate();
                         }
                         if(param1==1100)
                         {
                           _minHits=atoi(_minHitField->GetText());                         
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
                         if(param1==60)
                         {
                           _mainPad->cd();
                           TVirtualViewer3D *v;
                           v=_mainPad->GetViewer3D("ogl");
                         }
                         break;
   case kCM_CHECKBUTTON: if(param1==31)
                         {
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
                           if(_supportStructureButton->GetState()==kButtonDown)
                           {
                             _dataInterface->makeSupportStructuresVisible(true);
                           }
                           else
                           {
                             _dataInterface->makeSupportStructuresVisible(false);
                           }
                           drawEverything();
                         }
                         if(param1==33)
                         {
                           if(_outsideTracksButton->GetState()==kButtonDown)
                           {
                             double minx=_dataInterface->getTracksBoundary().minx;
                             double miny=_dataInterface->getTracksBoundary().miny;
                             double minz=_dataInterface->getTracksBoundary().minz;
                             double maxx=_dataInterface->getTracksBoundary().maxx;
                             double maxy=_dataInterface->getTracksBoundary().maxy;
                             double maxz=_dataInterface->getTracksBoundary().maxz;
                             _mainPad->GetView()->SetRange(minx,miny,minz,maxx,maxy,maxz);
                           }
                           else
                           {
                             double minx=_dataInterface->getHitsBoundary().minx;
                             double miny=_dataInterface->getHitsBoundary().miny;
                             double minz=_dataInterface->getHitsBoundary().minz;
                             double maxx=_dataInterface->getHitsBoundary().maxx;
                             double maxy=_dataInterface->getHitsBoundary().maxy;
                             double maxz=_dataInterface->getHitsBoundary().maxz;
                             _mainPad->GetView()->SetRange(minx,miny,minz,maxx,maxy,maxz);
                           }
                           _mainPad->GetView()->AdjustScales();
                           _mainPad->Modified();
                           _mainPad->Update();
                         }
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

  //resets all shapes
  std::list<boost::shared_ptr<VirtualShape> >::const_iterator iter;
  const std::list<boost::shared_ptr<VirtualShape> > &components=_dataInterface->getComponents();
  for(iter=components.begin(); iter!=components.end(); iter++)
  {
    (*iter)->start();
  }
  _mainPad->Modified();
  _mainPad->Update();
 
  if(_outsideTracksButton->GetState()==kButtonDown)
  {
    _timeStart=_dataInterface->getTracksBoundary().mint; 
    _timeStop=_dataInterface->getTracksBoundary().maxt;
  }
  else
  {
    _timeStart=_dataInterface->getHitsBoundary().mint; 
    _timeStop=_dataInterface->getHitsBoundary().maxt; 
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
    if(_saveAnim) combineAnimFiles();
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
  std::list<boost::shared_ptr<VirtualShape> >::const_iterator iter;
  const std::list<boost::shared_ptr<VirtualShape> > &components=_dataInterface->getComponents();
  for(iter=components.begin(); iter!=components.end(); iter++)
  {
    (*iter)->update(_timeCurrent);
  }
  
  if(!_clock)
  {
    char timeText[50];
    sprintf(timeText,"%+.4e ns",_timeCurrent);
    _clock = new TText(0.3,-0.9,timeText);
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
  std::list<boost::shared_ptr<VirtualShape> >::const_iterator iter;
  const std::list<boost::shared_ptr<VirtualShape> > &components=_dataInterface->getComponents();
  for(iter=components.begin(); iter!=components.end(); iter++)
  {
    (*iter)->update(NAN);
  }
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
