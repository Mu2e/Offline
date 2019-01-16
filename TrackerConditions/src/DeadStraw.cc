//
//
//

// C++ includes
#include <iostream>

// Mu2e includes
#include "TrackerConditions/inc/DeadStraw.hh"
#include "DataProducts/inc/StrawId.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  bool DeadStraw::isDead ( StrawId id,double hitpos) const {
    bool retval(false);
    auto ifnd = _deadstraws.find(Range(id));
    if(ifnd != _deadstraws.end()){
      if(ifnd->_range > 0)
        retval = fabs(hitpos) < ifnd->_range;
      else
        retval = fabs(hitpos) > -ifnd->_range;
    }
    return retval;
  }

  void DeadStraw::print( ostream& out) const{
    for( auto idead: _deadstraws) {
      out << "Straw " << idead._strawId 
          << " is dead for distances < " << idead._range << endl;
    }
  }
  
} // namespace mu2e
