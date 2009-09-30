/*------------------------------------------------------------

  A filter module that keeps only selected events.  It always 
  keeps event 3.  It never keeps event 4.  Otherwise it 
  selects events with either even or odd event numbers; this 
  choice is controlled by a run time configuration parameter.

  $Id: Ex02SelectEvents_plugin.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
  $Author: kutschke $
  $Date: 2009/09/30 22:57:47 $
   
  Original author Rob Kutschke


  -----------------------------------------------------------*/

// C++ includes
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include <boost/shared_ptr.hpp>


namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class Ex02SelectEvents : public edm::EDFilter {
  public:
    explicit Ex02SelectEvents(edm::ParameterSet const& pset):
      _keepOddOrEven(pset.getUntrackedParameter<int>("keepOddOrEven",1)){
    }
    virtual ~Ex02SelectEvents() { }
    virtual bool filter(edm::Event& e, edm::EventSetup const& c);
    
  private:

    // Control parameter: 1 to select odd numbered events; 
    //                    else select even numbered event.
    int _keepOddOrEven;

  };

  bool Ex02SelectEvents::filter(edm::Event& e, edm::EventSetup const&) {

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
DEFINE_FWK_MODULE(Ex02SelectEvents);
