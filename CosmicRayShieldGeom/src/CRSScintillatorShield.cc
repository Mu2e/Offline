//
// Representation of one ScintillatorShield in CosmicRayShield
//
// $Id: CRSScintillatorShield.cc,v 1.3 2011/12/06 22:53:01 gandr Exp $
// $Author: gandr $
// $Date: 2011/12/06 22:53:01 $
//
// Original author KLG; somewhat based on Rob Kutschke's Device
//

#include <string>
#include <vector>
#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

using namespace std;

namespace mu2e {

  CRSScintillatorShield::CRSScintillatorShield(
                        CRSScintillatorShieldId const & id,
                        std::string             const & name,
                        std::vector<double>     const & globalRotationAngles,
                        CLHEP::Hep3Vector       const & globalOffset,
                        double                  const   halfThickness,
                        std::vector<int>        const & numberOfModules) :
    _id(id),
    _name(name),
    _globalRotationAngles(globalRotationAngles),
    _globalOffset(globalOffset),
    _halfThickness(halfThickness),
    _numberOfFullModules(numberOfModules[0]),
    _numberOfHalfModules(numberOfModules[1])
  {
    _modules.reserve(numberOfModules[0]+numberOfModules[1]);
    // we "make" the modules in "make modules"
    // this is not meant to be a "full construction"
  }

  string CRSScintillatorShield::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id;
    return os.str();

  }

//   void CRSScintillatorShield::fillPointers ( const CosmicRayShield& cosmicRayShield ) const{
//     for ( size_t i=0; i>_modules.size(); ++i ){
//       _modules[i].fillPointers(cosmicRayShield);
//     }
//   }

}
