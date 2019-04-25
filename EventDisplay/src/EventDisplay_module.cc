//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: EventDisplay_module.cc,v 1.24 2014/04/18 16:46:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/04/18 16:46:25 $
//

#include <iostream>
#include <memory>
#include <string>

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "TApplication.h"
#include "TGMsgBox.h"

#include "EventDisplay/src/EventDisplayFrame.h"
#include "EventDisplay/src/RootFileManager.h"

using namespace CLHEP;
#include "RecoDataProducts/inc/KalRepCollection.hh"

namespace mu2e
{
  class EventDisplay : public art::EDAnalyzer
  {
    public:
    explicit EventDisplay(fhicl::ParameterSet const &pset);
    virtual ~EventDisplay() { }
    virtual void beginJob();
    void endJob();
    void analyze(const art::Event& e);

    private:
    template<class collectionType> void checkMinimumHits(const art::Event &event,
                                                         const std::string &classNameToCheck,
                                                         const mu2e_eventdisplay::EventDisplayFrame *_frame,
                                                         bool &showEvent);
    void checkMinimumHitsKalman(const art::Event &event, 
                                const mu2e_eventdisplay::EventDisplayFrame *_frame, 
                                bool &showEvent);
    mu2e_eventdisplay::EventDisplayFrame *_frame;
    bool _firstLoop;

    fhicl::ParameterSet _pset;
  };

  EventDisplay::EventDisplay(fhicl::ParameterSet const &pset)
    :
    art::EDAnalyzer(pset),
    _frame(0),
    _firstLoop(true),
    _pset(pset)
  {}

  void EventDisplay::beginJob()
  {
    //Bare pointer are needed for gApplication and _mainframe to avoid
    //that the destructor for these two objects gets called, i.e.
    //there cannot be a "delete" for these two objects.
    //The reason is that at least one of the ROOT destructors
    //(which a destructor call for these two objects would trigger)
    //has a bug (perhaps deletes things twice), so that the program
    //would crash.
    if (!gApplication) gApplication = new TApplication("EventDisplay",0,0);

    //don't create the eventdisplay here to avoid an empty (and not-responding)
    //eventdisplay window, because it usually takes a while until the first event gets pushed through
  }

  void EventDisplay::analyze(const art::Event& event)
  {
    TVirtualPad *temp_pad=gPad;
    TDirectory  *temp_dir=gDirectory;
    if(_firstLoop)
    {
      int x,y;
      unsigned int width,height;
      gVirtualX->GetWindowSize(gClient->GetRoot()->GetId(),x,y,width,height);
      width-=30;
      height-=70;
      _frame = new mu2e_eventdisplay::EventDisplayFrame(gClient->GetRoot(), width, height, _pset);
      if(!_frame->isClosed()) _frame->fillGeometry();
    }
    if(!_frame->isClosed())
    {
      bool findEvent=false;
      int eventToFind=_frame->getEventToFind(findEvent);
      if(findEvent)
      {
        int eventNumber=event.id().event();
        if(eventNumber==eventToFind) _frame->setEvent(event);
        else std::cout<<"event skipped, since this is not the event we are looking for"<<std::endl;
      }
      else
      {
        bool showEvent=true;
        checkMinimumHits<mu2e::StepPointMCCollection>(event, "<mu2e::StepPointMC>", _frame, showEvent);
        checkMinimumHits<mu2e::StepPointMCCollection>(event, "<StepPointMC>", _frame, showEvent);
        checkMinimumHits<mu2e::StrawHitCollection>(event, "<mu2e::StrawHit>", _frame, showEvent);
        checkMinimumHits<mu2e::StrawHitCollection>(event, "<StrawHit>", _frame, showEvent);
        checkMinimumHitsKalman(event, _frame, showEvent);
        if(showEvent) _frame->setEvent(event,_firstLoop);
      }
    }

    _firstLoop=false;
    if(temp_pad) temp_pad->cd(); else gPad=nullptr;
    if(temp_dir) temp_dir->cd(); else gDirectory=nullptr;

    if(_frame->isClosed()) 
    {
      throw cet::exception("CONTROL")<<"QUIT\n";
    }
  }

  template<class collectionType>
  void EventDisplay::checkMinimumHits(const art::Event &event, const std::string &classNameToCheck,
                                      const mu2e_eventdisplay::EventDisplayFrame *_frame, bool &showEvent)
  {
    std::string className, moduleLabel, productInstanceName;
    bool hasSelectedHits=_frame->getSelectedHitsName(className, moduleLabel, productInstanceName);
    if(hasSelectedHits && className.find(classNameToCheck)!=std::string::npos)
    {
      art::Handle<collectionType> hits;
      if(event.getByLabel(moduleLabel,productInstanceName,hits))
      {
        if(static_cast<int>(hits->size()) < _frame->getMinimumHits())
        {
          std::cout<<"event skipped, since it doesn't have enough hits"<<std::endl;
          showEvent=false;
        }
      }
    }
  }

  void EventDisplay::checkMinimumHitsKalman(const art::Event &event,
                                            const mu2e_eventdisplay::EventDisplayFrame *_frame, bool &showEvent)
  {
    std::string className, moduleLabel, productInstanceName;
    bool hasSelectedHits=_frame->getSelectedHitsName(className, moduleLabel, productInstanceName);
    if(hasSelectedHits && ((className.find("<KalRep>")!=std::string::npos) || (className.find("<mu2e::KalRep>")!=std::string::npos)))
    {
      art::Handle<mu2e::KalRepCollection> kalmantrackCollection;
      if(event.getByLabel(moduleLabel,productInstanceName,kalmantrackCollection))
      {
        int numberHits=0;
        for(unsigned int i=0; i<kalmantrackCollection->size(); i++)
        {
          KalRep const* particle = kalmantrackCollection->get(i);
          TrkHitVector const& hots = particle->hitVector();
	  numberHits+=hots.size();
        }
        if(numberHits < _frame->getMinimumHits())
        {
          std::cout<<"event skipped, since it doesn't have enough hits"<<std::endl;
          showEvent=false;
        }
      }
    }
  }

  void EventDisplay::endJob()
  {
    _frame->getRootFileManager()->write();
    if(!_frame->isClosed())
    {
      bool findEvent=false;
      int eventToFind=_frame->getEventToFind(findEvent);
      if(findEvent)
      {
        char msg[300];
        sprintf(msg,"The end of file has been reached, but the event #%i has not been found.",eventToFind);
        new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),"Event Not Found",msg,kMBIconExclamation,kMBOk);
      }
      _frame->CloseWindow();
    }
  }
}

using mu2e::EventDisplay;
DEFINE_ART_MODULE(EventDisplay);
