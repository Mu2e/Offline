//
// An anti-leak system to aid in using G4 from the Mu2e framework.
//
//
// Original author Rob Kutschke
//
// See header file for details.
//

#include <iomanip>

#include "G4Helper/inc/AntiLeakRegistry.hh"

namespace mu2e {

  void AntiLeakRegistry::clear(){
    _mapOfLists.clear();
    _maxNameLength = 0;
  }


  // Print information about types and numbers of saved objects.
  void AntiLeakRegistry::print( std::ostream& ost ) const{

    ost << "\nMu2e Registry for G4 objects " << std::endl;
    ost << "Number of entries: " << _mapOfLists.size() << std::endl;

    // I hope we never have more than 999,999 entries in any map.
    static const size_t sizeFieldWidth(6);

    for ( ListMap::const_iterator i = _mapOfLists.begin();
          i != _mapOfLists.end(); ++i ){
      ost << "  Data type: "
          << std::setw(_maxNameLength)
          << i->first
          << " Number of entries: "
          << std::setw(sizeFieldWidth)
          << i->second->size()
          << std::endl;
    }

  } // end AntiLeakRegistry::print

} // end namespace mu2e

