//
// Representation of one CRSScintillatorBar in CosmicRayShield
//
// $Id: CRSScintillatorBar.cc,v 1.1 2011/03/09 19:41:20 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:41:20 $
//
// Original author KLG
//

#include <sstream>
#include "CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e {

  CRSScintillatorBar::CRSScintillatorBar():
    _id(CRSScintillatorBarId()),
    _index(CRSScintillatorBarIndex(0)),
    _localOffset (CLHEP::Hep3Vector(0.,0.,0.)),
    _globalRotationAngles(std::vector<double>(2)),
    _globalOffset(CLHEP::Hep3Vector(0.,0.,0.))
  {}

  CRSScintillatorBar::CRSScintillatorBar(CRSScintillatorBarId const& id,
                                         CRSScintillatorBarIndex const& index
                                         ):
    _id(id),
    _index(index),
    _localOffset (CLHEP::Hep3Vector(0.,0.,0.)),
    _globalRotationAngles(std::vector<double>(2)),
    _globalOffset(CLHEP::Hep3Vector(0.,0.,0.))
  {}

  // Constructor using bar normal unit vector.
  CRSScintillatorBar::CRSScintillatorBar(
                                         CRSScintillatorBarId const& id,
                                         CRSScintillatorBarIndex const& index,
                                         CLHEP::Hep3Vector const& localOffset,
                                         std::vector<double> const & globalRotationAngles,
                                         CLHEP::Hep3Vector const& globalOffset
                                         ):
    _id(id),
    _index(index),
    _localOffset(localOffset),
    _globalRotationAngles(globalRotationAngles),
    _globalOffset(globalOffset)
  {
    // we need to get or calculate the global offset
  }

//   void CRSScintillatorBar::fillPointers ( const CosmicRayShield& cosmicRayShield ) const{
//     _detail = &cosmicRayShield.getCRSScintillatorBarDetails().at(_detailIndex);
//   }

//   bool CRSScintillatorBar::isNearestNeighbour( CRSScintillatorBarIndex idx ) const{

//     for ( vector<CRSScintillatorBarIndex>::const_iterator i=_nearestByIndex.begin(),
//             e=_nearestByIndex.end(); 
//           i<e; ++i ){
//       if ( *i == idx ) return true;
//     }

//     return false;
//   }

  std::string CRSScintillatorBar::name( std::string const& base ) const{
    std::ostringstream os;

    os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber() << "_"
       << _id.getLayerNumber()  << "_"
       << _id.getBarNumber();
    return os.str();

  }

} // namespace mu2e

