//
// An enum-matched-to-names class used to indicate the level of compression
// requested. The class contains enum entries for the different possible
// levels of compression.
//
//
// Original author Rob Kutschke
// Modified for Compression: Andy Edmonds, Dec 2020
//

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Compression/inc/CompressionLevel.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {

  const char* CompressionLevel::_name[] = { COMPRESSIONLEVEL_NAMES };

  BOOST_STATIC_ASSERT(sizeof(CompressionLevel::_name)/sizeof(char*) == CompressionLevel::lastEnum);

  void CompressionLevel::printAll( std::ostream& ost){
    ost << "List of allowed compression levels: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(3) << i << " " << _name[i] << std::endl;
    }
  }

  CompressionLevel CompressionLevel::findByName ( std::string const& name){

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<size(); ++i ){
      if ( _name[i] == name ){
        return CompressionLevel(enum_type(i));
      }
    }
    return CompressionLevel(unknown);
  }
}
