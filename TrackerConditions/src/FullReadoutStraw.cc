//
//
//


#include <iostream>
#include "TrackerConditions/inc/FullReadoutStraw.hh"
#include "DataProducts/inc/StrawId.hh"
#include <iostream>

using namespace std;

namespace mu2e {

  void FullReadoutStraw::print( ostream& out) const{
    for( auto is: _straws) {
      out << "Straw " << is << " has full readout" << endl;
    }
  }
  
} // namespace mu2e
