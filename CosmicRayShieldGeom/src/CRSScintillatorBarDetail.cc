//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.cc,v 1.4 2012/03/29 22:56:21 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/03/29 22:56:21 $
//
// Original author KLG; somewhat based on Rob Kutschke StrawDetail
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"

using namespace std;

namespace mu2e {

  CRSScintillatorBarDetail::CRSScintillatorBarDetail( int const id,
                                                      std::vector<std::string> const& materialNames,
                                                      std::vector<double> const& halfLengths
                                                      ) :
    _id(id),
    _materialNames(materialNames),
    _halfLengths(halfLengths)
  {}


  // Construct a string containing the bar Id.
  std::string CRSScintillatorBarDetail::name( std::string const& base ) const{
    std::ostringstream os;
    os << base
       << _id;
    return os.str();
  }


} // namespace mu2e
