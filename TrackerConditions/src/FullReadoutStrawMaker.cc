//
//
//

#include "TrackerConditions/inc/FullReadoutStraw.hh"
#include "TrackerConditions/inc/FullReadoutStrawConfig.hh"
#include "TrackerConditions/inc/FullReadoutStrawMaker.hh"
#include "DataProducts/inc/StrawId.hh"
#include <iostream>
#include <sstream>
#include <memory>

using namespace std;

namespace mu2e {


  FullReadoutStraw::ptr_t FullReadoutStrawMaker::fromFcl() {

    // the list to be filled
    FullReadoutStraw::vec_t straws;

    for(auto ss : _config.straws()) straws.emplace_back(StrawId(ss));

    // make shared_ptr to the FullReadoutStraw object on the heap
    FullReadoutStraw::ptr_t ptr = std::make_shared<FullReadoutStraw>(straws);

    if ( _config.verbose() > 0 ) {
      cout << "FullReadoutStrawMaker created a list of " << straws.size() 
	   << " straws" << endl;
    }
    
    return ptr;
  }

  FullReadoutStraw::ptr_t FullReadoutStrawMaker::fromDb() {
    return fromFcl();
  }


} // namespace mu2e
