//
// An enum-matched-to-names class for STM channels
// based on GenId.cc
#include <iomanip>
#include <iostream>

#include "Offline/DataProducts/inc/STMChannel.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {

  const char* STMChannel::_name[] = { STMCHANNEL_NAMES };

  BOOST_STATIC_ASSERT(sizeof(STMChannel::_name)/sizeof(char*) == STMChannel::lastEnum);

  void STMChannel::printAll( std::ostream& ost){
    ost << "List of STM channels: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(2) << i << " " << _name[i] << std::endl;
    }
  }

  STMChannel STMChannel::findByName ( std::string const& name){

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<size(); ++i ){
      if ( _name[i] == name ){
        return STMChannel(enum_type(i));
      }
    }
    return STMChannel(unknown);
  }


}
