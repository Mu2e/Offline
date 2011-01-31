//
// Module which starts the event display, and transmits the data of each event to the event display.
//
// $Id: EventDisplay_plugin.cc,v 1.2 2011/01/31 01:03:53 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/01/31 01:03:53 $
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

#include "TApplication.h"

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
    _frame = new mu2e_eventdisplay::EventDisplayFrame(gClient->GetRoot(), 800, 550);
  }

  void EventDisplay::analyze(const edm::Event& event, edm::EventSetup const&) 
  {
    if(_firstLoop)
    {
      _firstLoop=false;
      if(!_frame->isClosed()) _frame->fillGeometry();
    }
    if(!_frame->isClosed())
    {
/*This section needs to be completely rewritten */
      edm::Handle<mu2e::StepPointMCCollection> hits;
      std::string _g4ModuleLabel = "g4run";
      std::string _trackerStepPoints = "tracker"; //TODO: this may not always be correct 
                                         //in the future: let user decide via display,
                                         //and get this info with a get function from _frame
                                         //however, something (or perhaps everything) needs
                                         //to be selected by default, otherwise all events 
                                         //will be skipped, and the user will never have a 
                                         //chance to select anything
      if(event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits))  
      {
        int numberHits=hits->size();
        if(numberHits > _frame->getMinimumHits()) _frame->fillEvent(event);
        else std::cout<<"event skipped, since it doesn't have enough hits"<<std::endl;
      }
/*-----------------------------------*/
    }
  }

  void EventDisplay::endJob()
  {
    if(!_frame->isClosed()) _frame->CloseWindow();
  }
}

using mu2e::EventDisplay;
DEFINE_FWK_MODULE(EventDisplay);
