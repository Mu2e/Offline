//
// An enum-matched-to-names class for Virtual Detector Id's.
//
// Original author Rob Kutschke

#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "cetlib_except/exception.h"

#include <boost/static_assert.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

  VirtualDetectorId::VirtualDetectorId( std::string const& name,
                                        bool throwIfUnknown,
                                        bool throwIfUndefined)
    :_id(unknown)
  {

    auto b = names().begin();
    auto e = names().end();
    auto i = find( b, e , name );

    // Found it.
    if ( i != e ){
      _id = static_cast<enum_type>( i - b );
      if ( _id == unknown && throwIfUnknown ) {
        throw cet::exception("BADCONFIG") << "VirtualDetectorId::unknown is not allowed at this time";
      }
      return;
    }

    // Did not find it.
    if ( throwIfUndefined ) {
      throw cet::exception("BADCONFIG") << "VirtualDetectorID: invalid enum name : " << name;
    }

  }


  void VirtualDetectorId::printAll( std::ostream& ost){
    ost << "List of Virtual Detector Id codes: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      VirtualDetectorId id(i);
      ost << setw(3) << i << " " << name(id.id()) << std::endl;
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
      throw cet::exception("OutOfRange") << "Invalid VirtualDetectorID::enum_type : " << id;
    }
    return true;
  }

  std::ostream& operator<<(std::ostream& ost,
                           const VirtualDetectorId& id ){
    ost << "( "
        << id.id() << ": "
        << id.name()
        << " )";
    return ost;
  }



}
