//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.cc,v 1.3 2011/05/18 15:47:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:47:40 $
//
// Original author KLG; somewhat based on Rob Kutschke StrawDetail
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarIndex.hh"

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
