//
// construct the TrackerStatus object from fcl configuration
//

// C++ includes
#include <iostream>
#include <sstream>
#include <memory>

// Mu2e includes
#include "TrackerConditions/inc/TrackerStatusMaker.hh"

using namespace std;

namespace mu2e {

  TrackerStatus::ptr_t TrackerStatusMaker::fromFcl() {

    TrackerStatus::estat_t estats;
    for(auto const& estat : config_.status()) {
      estats.emplace(StrawId(std::get<0>(estat)),
	  StrawIdMask(std::get<1>(estat)), 
	  StrawStatus(std::get<2>(estat)));
    }
    auto const& settings = config_.settings();
    if ( settings.verbose() > 0 ) {
      cout << "TrackerStatus created with " << estats.size() 
	<< " elements " << endl;
      if(settings.verbose() > 1){
	for(auto const& estat : estats) {
	  std::cout << "Tracker Element Id " << estat.sid_
	    << " level " << estat.mask_.levelName()  
	    << " mask " << estat.mask_.mask()  
	    << " status " << estat.status_ << std::endl;
	}
      }
    }
    // make shared_ptr to the TrackerStatus object on the heap
    return std::make_shared<TrackerStatus>(estats);
  }

  TrackerStatus::ptr_t TrackerStatusMaker::fromDb() {
    return fromFcl(); // TODO!!!
  }

  } // namespace mu2e
