//
// An enum-matched-to-names class for magnetic field types.
//
//
// $Id: BFMapType.cc,v 1.2 2011/05/18 02:27:14 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:14 $
//
// Original author Rob Kutschke

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "BFieldGeom/inc/BFMapType.hh"

#include <boost/static_assert.hpp>

using namespace std;

namespace mu2e {

  BFMapType::BFMapType( int id):
    _id(static_cast<enum_type>(id)){
    if ( !isValid(_id) ){
      ostringstream os;
      os << "BFMapType: identifier is out of range: " << id
         << "  Allowed range: [" << unknown << "," << lastEnum << ")";
      throw std::out_of_range(os.str());
    }
  }

  std::string const& BFMapType::name( enum_type id){

    static std::vector<string> names;

    // Intiailize the vector on the first call.
    if ( names.size() == 0 ){
      const char* _xname[] = { BFMAPTYPE_NAMES };
      const size_t n(sizeof(_xname)/sizeof(char*));
      BOOST_STATIC_ASSERT(n == BFMapType::lastEnum);
      for ( size_t i = 0; i<n; ++i){
        names.push_back(_xname[i]);
      }
    }

    // Belt and suspenders.
    return names.at(id);
  }

  void BFMapType::printAll( std::ostream& ost){
    ost << "List of BFMapType id codes: " << endl;
    for ( int i=unknown; i<lastEnum; ++i){
      enum_type j = static_cast<enum_type>(i);
      ost << setw(2) << i << " " << name(j) << std::endl;
    }
  }

} // end namespace mu2e
