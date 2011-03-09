//
// Representation of common properties of the Scintillator Bars etc...
//
// $Id: CRSScintillatorBarDetail.cc,v 1.1 2011/03/09 19:41:31 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:41:31 $
//
// Original author KLG; somewhat based on Rob Kutschke StrawDetail
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorBarDetail.hh"
#include "CosmicRayShieldGeom/inc/CRSScintillatorBarIndex.hh"

using namespace std;

namespace mu2e {

  CRSScintillatorBarDetail::CRSScintillatorBarDetail( int32_t const id,
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
