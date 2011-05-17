/*------------------------------------------------------------

  A filter module that keeps only selected events.  It always 
  keeps event 3.  It never keeps event 4.  Otherwise it 
  selects events with either even or odd event numbers; this 
  choice is controlled by a run time configuration parameter.

  $Id: Ex02SelectEvents_plugin.cc,v 1.2 2011/05/17 15:36:00 greenc Exp $
  $Author: greenc $
  $Date: 2011/05/17 15:36:00 $
   
  Original author Rob Kutschke


  -----------------------------------------------------------*/

// C++ includes
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include <boost/shared_ptr.hpp>


namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class Ex02SelectEvents : public art::EDFilter {
  public:
    explicit Ex02SelectEvents(fhicl::ParameterSet const& pset):
      _keepOddOrEven(pset.get<int>("keepOddOrEven",1)){
    }
    virtual ~Ex02SelectEvents() { }
    virtual bool filter(art::Event& e, art::EventSetup const& c);
    
  private:

    // Control parameter: 1 to select odd numbered events; 
    //                    else select even numbered event.
    int _keepOddOrEven;

  };

  bool Ex02SelectEvents::filter(art::Event& e, art::EventSetup const&) {

    // EventSetup is a cms leftover that we do not use.
    
    // Get event number from the event.
    int event = e.id().event();
    
    // Always keep event 3.
    if ( event == 3 ) return true;

    // Always discard event 4.
    if ( event == 4 ) return false;

    // Is this an odd numbered event?
    bool isOdd = ((event % 2) != 0);
  
    // Keep only events with odd or even event numbers, as 
    // controled by the parameter from the ParameterSet.
    if ( _keepOddOrEven == 1){
      return isOdd;
    } else{
      return !isOdd;
    }
  
  }

}


using mu2e::Ex02SelectEvents;
DEFINE_ART_MODULE(Ex02SelectEvents);
