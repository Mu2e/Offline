//
//
//

// C++ includes
#include <iostream>
#include <sstream>
#include <memory>

// Mu2e includes
#include "TrackerConditions/inc/DeadStraw.hh"
#include "TrackerConditions/inc/DeadStrawMaker.hh"
#include "DataProducts/inc/StrawId.hh"

using namespace std;

namespace mu2e {


  DeadStraw::ptr_t DeadStrawMaker::fromFcl() {

    // the list to be filled
    DeadStraw::set_t deadstraws;

    for(auto ii : _config.deadPlanes()) addDeadPlane(PlaneId(ii,0,0),deadstraws);
    for(auto ss : _config.deadPanels()) addDeadPanel(PanelId(ss),deadstraws);
    for(auto ss : _config.deadStraws()) addDeadStraw(StrawId(ss),deadstraws);
    for(auto ss : _config.partlyDeadStraws()) {
      // interpret "StrawID range" strings like "11_2_1 10.0"
      istringstream dstrings(ss);
      double range(-1.0);
      string sidname;
      dstrings >> sidname >> range;
      // check
      if(range < 0.0)
	throw cet::exception("DEAD_STRAW_FCL_RANGE")
	  << "DeadStraw: expected StrawId and range but got " << ss << endl;
      addDeadStraw(StrawId(sidname),deadstraws,range);
    }

    // make shared_ptr to the DeadStraw object on the heap
    DeadStraw::ptr_t ptr = std::make_shared<DeadStraw>(deadstraws);

    if ( _config.verbose() > 0 ) {
      cout << "DeadStrawMaker created " << deadstraws.size() 
	   << " dead straws" << endl;
    }

    
    return ptr;
  }

  DeadStraw::ptr_t DeadStrawMaker::fromDb() {
    return fromFcl();
  }

  void DeadStrawMaker::addDeadPlane( PlaneId const& id, DeadStraw::set_t& dead){
    for(uint16_t panel=0; panel<StrawId::_npanels; panel++) 
      addDeadPanel(PanelId(id.getPlane(),panel,0),dead);
  }

  void DeadStrawMaker::addDeadPanel( PanelId const& id, DeadStraw::set_t& dead){
    for(uint16_t straw=0; straw<StrawId::_nstraws; straw++) 
      addDeadStraw(StrawId(id.getPlane(),id.getPanel(),straw),dead);
  }
 
  void DeadStrawMaker::addDeadStraw( StrawId const& id, DeadStraw::set_t& dead, 
				     double range){
    if(range>0.0) {
      dead.emplace(id,range);
    } else {
      dead.emplace(id);
    }
  }


} // namespace mu2e
