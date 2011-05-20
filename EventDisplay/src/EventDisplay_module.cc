//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: EventDisplay_module.cc,v 1.5 2011/05/20 20:18:23 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:18:23 $
//

#include <iostream>
#include <memory>
#include <string>

#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Persistency/Common/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "TApplication.h"
#include "TGMsgBox.h"

#include "EventDisplayFrame.h"

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
    mu2e_eventdisplay::EventDisplayFrame *_frame;
    bool _firstLoop;

    fhicl::ParameterSet const &_pset;
  };

  EventDisplay::EventDisplay(fhicl::ParameterSet const &pset)
    :
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
    if(_firstLoop)
    {
      _frame = new mu2e_eventdisplay::EventDisplayFrame(gClient->GetRoot(), 800, 550, _pset);
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
        checkMinimumHits<mu2e::StepPointMCCollection>(event, "std::vector<mu2e::StepPointMC>", _frame, showEvent);
        checkMinimumHits<mu2e::StrawHitCollection>(event, "std::vector<mu2e::StrawHit>", _frame, showEvent);
        if(showEvent) _frame->setEvent(event,_firstLoop);
      }
    }
    _firstLoop=false;
  }

  template<class collectionType>
  void EventDisplay::checkMinimumHits(const art::Event &event, const std::string &classNameToCheck,
                                      const mu2e_eventdisplay::EventDisplayFrame *_frame, bool &showEvent)
  {
    std::string className, moduleLabel, productInstanceName;
    bool hasSelectedHits=_frame->getSelectedHitsName(className, moduleLabel, productInstanceName);
    if(hasSelectedHits && className.compare(classNameToCheck)==0)
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

  void EventDisplay::endJob()
  {
    if(!_frame->isClosed())
    {
      bool findEvent=false;
      int eventToFind=_frame->getEventToFind(findEvent);
      if(findEvent)
      {
        char msg[300];
        sprintf(msg,"The end of file has been reached, but the event #%i has not been found.",eventToFind);
        TGMsgBox *eventNotFoundBox;
        eventNotFoundBox = new TGMsgBox(gClient->GetRoot(),gClient->GetRoot(),"Event Not Found",msg,kMBIconExclamation,kMBOk);
      }
      _frame->CloseWindow();
    }
  }
}

using mu2e::EventDisplay;
DEFINE_ART_MODULE(EventDisplay);
