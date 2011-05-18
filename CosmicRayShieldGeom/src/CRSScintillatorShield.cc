//
// Representation of one ScintillatorShield in CosmicRayShield
//
// $Id: CRSScintillatorShield.cc,v 1.2 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
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
                        CLHEP::Hep3Vector       const & localOffset,
                        std::vector<double>     const & globalRotationAngles,
                        CLHEP::Hep3Vector       const & globalOffset,
                        double                  const   halfThickness,
                        std::vector<int>        const & numberOfModules) :
    _id(id),
    _name(name),
    _localOffset(localOffset),
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
