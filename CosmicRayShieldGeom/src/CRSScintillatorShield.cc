//
// Representation of one ScintillatorShield in CosmicRayShield
//
// $Id: CRSScintillatorShield.cc,v 1.4 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG; somewhat based on Rob Kutschke's Device
//

#include <string>
#include <vector>
#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

using namespace std;

namespace mu2e {

  CRSScintillatorShield::CRSScintillatorShield(CRSScintillatorShieldId const & id,
                                               std::string const & name,
                                               CRSScintillatorBarDetail const &barDetails) : 
    _id(id),
    _name(name),
    _barDetails(barDetails)
  {
  }

  string CRSScintillatorShield::name( string const& base ) const
  {
    ostringstream os;

    os << base
       << _id;
    return os.str();
  }

}
