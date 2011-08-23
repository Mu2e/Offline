//
// An enum-matched-to-names class for Virtual Detector Id's.
//
// $Id: VirtualDetectorId.cc,v 1.1 2011/08/23 21:48:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/08/23 21:48:43 $
//
// Original author Rob Kutschke

#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include <boost/static_assert.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  void VirtualDetectorId::printAll( std::ostream& ost){
    ost << "List of Virtual Detector Id codes: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      VirtualDetectorId id(i);
      ost << setw(2) << i << " " << name(id.id()) << std::endl;
    }
  }

  std::string const& VirtualDetectorId::name( enum_type id){
    return names().at(id);
  }

  std::vector<std::string> const& VirtualDetectorId::names(){

    static vector<string> nam;

    if ( nam.empty() ){
      const char* tmp[] = { VIRTUALDETECTORID_NAMES };
      BOOST_STATIC_ASSERT(sizeof(tmp)/sizeof(char*) == lastEnum);
      nam.insert( nam.begin(), tmp, tmp+lastEnum);
    }

    return nam;
  }

  bool VirtualDetectorId::isValidorThrow( enum_type id){
    if ( !isValid(id) ){
      ostringstream os;
      os << "Invalid VirtualDetectorId::enum_type : "
         << id;
      throw std::out_of_range( os.str() );
    }
    return true;
  }




}
