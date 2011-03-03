//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: EventDisplay_plugin.cc,v 1.6 2011/03/03 17:24:10 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/03/03 17:24:10 $
//

#include <iostream>
#include <string>
#include <memory>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"

#include "TApplication.h"
#include "TGMsgBox.h"

#include "EventDisplayFrame.h"

namespace mu2e 
{
  class EventDisplay : public edm::EDAnalyzer 
  {
    public:
    explicit EventDisplay(edm::ParameterSet const&);
    virtual ~EventDisplay() { }
    virtual void beginJob(edm::EventSetup const&);
    void endJob();
    void analyze(const edm::Event& e, edm::EventSetup const&);

    private:
    template<class collectionType> void checkMinimumHits(const edm::Event &event, 
                                                         const std::string &classNameToCheck,
                                                         const mu2e_eventdisplay::EventDisplayFrame *_frame, 
                                                         bool &showEvent);
    mu2e_eventdisplay::EventDisplayFrame *_frame;
    bool _firstLoop;
  };

  EventDisplay::EventDisplay(edm::ParameterSet const& ) 
  {
    _firstLoop=true;
  }

  void EventDisplay::beginJob(edm::EventSetup const& )
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

  void EventDisplay::analyze(const edm::Event& event, edm::EventSetup const&) 
  {
    if(_firstLoop)
    {
      _frame = new mu2e_eventdisplay::EventDisplayFrame(gClient->GetRoot(), 800, 550);
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
  void EventDisplay::checkMinimumHits(const edm::Event &event, const std::string &classNameToCheck,
                                      const mu2e_eventdisplay::EventDisplayFrame *_frame, bool &showEvent)
  {
    std::string className, moduleLabel, productInstanceName;
    bool hasSelectedHits=_frame->getSelectedHitsName(className, moduleLabel, productInstanceName);
    if(hasSelectedHits && className.compare(classNameToCheck)==0)
    {
      edm::Handle<collectionType> hits;
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
DEFINE_FWK_MODULE(EventDisplay);
