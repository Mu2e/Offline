//
// Representation of one Scintillator Module in  CosmicRayShield
//
//
// $Id: CRSScintillatorModule.cc,v 1.1 2011/03/09 19:42:24 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:42:24 $
//
// Original author KLG based on Rob Kutschke's Sector
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"

using namespace std;

namespace mu2e {

    CRSScintillatorModule::CRSScintillatorModule():
      _id(CRSScintillatorModuleId(-1,-1)),
      _nBarsPerLayer(0),
      _localOffset( CLHEP::Hep3Vector(0.,0.,0.)),
      _globalRotationAngles(std::vector<double>(2)),
      _globalOffset(CLHEP::Hep3Vector(0.,0.,0.))
    {};

    CRSScintillatorModule::CRSScintillatorModule( CRSScintillatorModuleId const& id, 
                                                  int const nBarsPerLayer ):
      _id(id),
      _nBarsPerLayer(nBarsPerLayer)
    {};

    CRSScintillatorModule::CRSScintillatorModule(
                                                 CRSScintillatorModuleId const & id,
                                                 int               const nBarsPerLayer,
                                                 CLHEP::Hep3Vector const & localOffset,
                                                 std::vector<double> const & globalRotationAngles,
                                                 CLHEP::Hep3Vector const & globalOffset // offset in World
                                                 ):
      _id(id),
      _nBarsPerLayer(nBarsPerLayer),
      _localOffset(localOffset),
      _globalRotationAngles(globalRotationAngles),
      _globalOffset(globalOffset)
    {};


  string CRSScintillatorModule::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber();
    return os.str();

  }

//   void CRSScintillatorModule::fillPointers ( const CosmicRayShield& cosmicRayShield ) const {
//     for( size_t i=0; i<_layers.size(); ++i ){
//       _layers[i].fillPointers(cosmicRayShield);
//     }
//   }

} // end namespace mu2e
