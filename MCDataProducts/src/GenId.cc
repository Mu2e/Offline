a//
// An enum-matched-to-names class for generator Id's.
//
// Add a line that says fixme in a file that will trigger clang tidy.
// and FIXME
// and Fixme
//
// Original author Rob Kutschke
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "Offline/MCDataProducts/inc/GenId.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {

  const char* GenId::_name[] = { GENID_NAMES };

  BOOST_STATIC_ASSERT(sizeof(GenId::_name)/sizeof(char*) == GenId::lastEnum);

  void GenId::printAll( std::ostream& ost){
    ost << "List of Generator Id codes: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(2) << i << " " << _name[i] << std::endl;
    }
  }

  GenId GenId::findByName ( std::string const& name, bool throwIfUnknown, bool throwIfUndefined ){

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<size(); ++i ){
      if ( _name[i] == name ){
        if ( throwIfUnknown && enum_type(i) == unknown ){
          throw std::out_of_range( "GenId::unknown is not allowed at this time" );
        }
        return GenId(enum_type(i));
      }
    }

    if ( throwIfUndefined ){
      std::ostringstream os;
      os << "GenId::findByName invalid enum name : " << name;
      throw std::out_of_range( os.str() );
    }
    return GenId(unknown);  
  }

}

/*

Notes from Marc Paterno.  Other options for this class.

 One option is codegen.  I think this is the best solution but we
 don't happen to have handy technology right now.

 Once I have the char**, I can initialize a vector as:
 vector<string> vname(n, names, names+n);

 Since the text methods are for debug/info, the name
 method can return a string by copy.

 I tried to do the following in the header but it made
 a compiler error:
 static const char* name[] = {};

 The problem with making this a template is that if we add a new
 enum value this makes a new type - which will confuse the edm.

 Marc thinks that a default c'tor that assigns an initial value
 formally breaks POD-ness even though sizeof(GenID) = sizeof(int).


*/
