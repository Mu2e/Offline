//
// An enum-matched-to-names class for physics processes from G4,
// plus mu2e defined reasons for stopping particles.
//
// $Id: StoppingCode.cc,v 1.2 2010/12/17 22:21:43 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/12/17 22:21:43 $
//
// Original author Rob Kutschke

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "ToyDP/inc/StoppingCode.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {
  
  const char* StoppingCode::_name[] = { STOPPINGCODE_NAMES };

  BOOST_STATIC_ASSERT(sizeof(StoppingCode::_name)/sizeof(char*) == StoppingCode::lastEnum);

  void StoppingCode::printAll( std::ostream& ost){
    ost << "List of stopping codes from G4: " << endl;
    for ( int i=0; i<lastEnum; ++i){
      ost << setw(2) << i << " " << _name[i] << std::endl;
    }
  }

  StoppingCode::StoppingCode( int id ):
    _id(static_cast<enum_type>(id)){
    if ( !isValid() ){
      ostringstream os;
      os << "Invalid StoppingCode::enum_type: " << id;
      throw std::logic_error ( os.str());
    }
  }
  
  StoppingCode StoppingCode::findByName ( std::string const& name){

    // Size must be at least 2 (for unknown and lastEnum).
    for ( size_t i=0; i<size(); ++i ){
      if ( _name[i] == name ){
        return StoppingCode(i);
      }
    }
    return StoppingCode(unknown);
  }

  // Return a vector of the codes that are mu2e specific.
  std::vector<StoppingCode> StoppingCode::mu2eCodes(){
    std::vector<StoppingCode> codes;
    for ( size_t i=0; i<size(); ++i ){
      string name = _name[i];
      if ( name.substr(0,4) == "mu2e" ){
        codes.push_back(StoppingCode(i));
      }
    }
    return codes;
  }


}
