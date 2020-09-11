//
// An enum-matched-to-names class for the legal operating modes of a filter that selects events
// based on the occupancy of StepPointMCs in various detectors.
//
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/StepFilterMode.hh"

#include <boost/static_assert.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  void StepFilterMode::printAll( std::ostream& ost){
    ost << "List of names of instances for StepPointMCCollections Id codes: "
        << endl;
    for ( int i=0; i<lastEnum; ++i){
      StepFilterMode id{enum_type(i)};
      ost << setw(2) << i << " " << name(id.id()) << std::endl;
    }
  }

  std::string const& StepFilterMode::name( enum_type id){
    return names().at(id);
  }

  std::vector<std::string> const& StepFilterMode::names(){

    static vector<string> nam;

    if ( nam.empty() ){
      const char* tmp[] = { STEPFILTERMODE_NAMES };
      BOOST_STATIC_ASSERT(sizeof(tmp)/sizeof(char*) == lastEnum);
      nam.insert( nam.begin(), tmp, tmp+lastEnum);
    }

    return nam;
  }

  bool StepFilterMode::isValidorThrow( enum_type id){
    if ( !isValid(id) ){
      ostringstream os;
      os << "Invalid StepFilterMode::enum_type : "
         << id;
      throw std::out_of_range( os.str() );
    }
    return true;
  }

  StepFilterMode StepFilterMode::findByName( std::string const& name, bool throwIfUndefined ){
    std::vector<string> const& nam(names());

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<nam.size(); ++i ){
      if ( nam[i] == name ){
        if ( i == unknown && throwIfUndefined ){
          throw std::out_of_range( "StepFilterMode the value \"unknown\" is not permitted in this context. " );
        }
        return StepFilterMode(enum_type(i));
      }
    }

    if ( throwIfUndefined ) {
      ostringstream os;
      os << "StepFilterMode unrecognized name: "
         << name;
      throw std::out_of_range( os.str() );
    }

    return StepFilterMode(unknown);
  } // end StepFilterMode::findByName

  // Get all values, as class instances; this includes the unknown value.
  std::vector<StepFilterMode> const& StepFilterMode::allValues(){

    static std::vector<StepFilterMode> values;
    for ( size_t i=unknown; i<lastEnum; ++i ){
      values.push_back(StepFilterMode(enum_type(i)));
    }
    return values;

  } // end StepFilterMode::allValues


} // end namespace mu2e
