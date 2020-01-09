//
//
// An enum-matched-to-names class for the names of the StepPointMC
// collections produced inside G4_module.
//
// $Id: StepInstanceName.cc,v 1.2 2012/06/09 20:04:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/09 20:04:51 $
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/StepInstanceName.hh"

#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  void StepInstanceName::printAll( std::ostream& ost){
    ost << "List of names of instances for StepPointMCCollections Id codes: "
        << endl;
    for ( int i=0; i<lastEnum; ++i){
      StepInstanceName id(i);
      ost << setw(2) << i << " " << name(id.id()) << std::endl;
    }
  }

  std::string const& StepInstanceName::name( enum_type id){
    return names().at(id);
  }

  std::vector<std::string> const& StepInstanceName::names() {
    static_assert(initializer_list{STEPINSTANCENAME_NAMES}.size() == lastEnum);
    static vector<string> const nam{STEPINSTANCENAME_NAMES};
    return nam;
  }

  bool StepInstanceName::isValidorThrow( enum_type id){
    if ( !isValid(id) ){
      ostringstream os;
      os << "Invalid StepInstanceName::enum_type : "
         << id;
      throw std::out_of_range( os.str() );
    }
    return true;
  }

  StepInstanceName StepInstanceName::findByName( std::string const& name, bool throwIfUndefined ){
    std::vector<string> const& nam(names());

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<nam.size(); ++i ){
      if ( nam[i] == name ){
        if ( i == unknown && throwIfUndefined ){
          throw std::out_of_range( "StepInstanceName the value \"unknown\" is not permitted in this context. " );
        }
        return StepInstanceName(i);
      }
    }

    if ( throwIfUndefined ) {
      ostringstream os;
      os << "StepInstanceName unrecognized name: "
         << name;
      throw std::out_of_range( os.str() );
    }

    return StepInstanceName(unknown);
  } // end StepInstanceName::findByName

  // Get all values, as class instances; this includes the unknown value.
  std::vector<StepInstanceName> const& StepInstanceName::allValues(){

    static std::vector<StepInstanceName> values;
    for ( size_t i=unknown; i<lastEnum; ++i ){
      values.push_back(StepInstanceName(i));
    }
    return values;

  } // end StepInstanceName::allValues


} // end namespace mu2e
