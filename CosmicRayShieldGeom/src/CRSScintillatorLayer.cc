//
// Representation of one Scintillator Layer in CosmicRayShield
//
//
// $Id: CRSScintillatorLayer.cc,v 1.1 2011/03/09 19:41:38 genser Exp $
// $Author: genser $
// $Date: 2011/03/09 19:41:38 $
//
// Original author KLG based on Rob Kutschke's Layer
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#ifndef __CINT__ 

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e {

  CRSScintillatorLayer::CRSScintillatorLayer():
    _id(CRSScintillatorLayerId()),
    _nBars(0),
    _localOffset( CLHEP::Hep3Vector(0.,0.,0.)),
    _globalRotationAngles(std::vector<double>(2)),
    _globalOffset(CLHEP::Hep3Vector(0.,0.,0.))
  {}

  CRSScintillatorLayer::CRSScintillatorLayer(
                                             CRSScintillatorLayerId const&   id,
                                             int const nBars
                                             ):
    _id(id),
    _nBars(nBars)
  {}

  CRSScintillatorLayer::CRSScintillatorLayer(
                                             CRSScintillatorLayerId const&   id,
                                             int const nBars,
                                             CLHEP::Hep3Vector const& localOffset, // wrt Shield
                                             std::vector<double> const & globalRotationAngles,
                                             CLHEP::Hep3Vector const& globalOffset // wrt World
                                             ):
    _id(id),
    _nBars(nBars),
    _localOffset(localOffset),
    _globalRotationAngles(globalRotationAngles),
    _globalOffset(globalOffset)    
  {}

//   CRSScintillatorLayer::CRSScintillatorLayer(CRSScintillatorLayerId const& id ):
//     _id(id){
//   }

  string CRSScintillatorLayer::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber() << "_"
       << _id.getLayerNumber();
    return os.str();

  }

//   void CRSScintillatorLayer::fillPointers ( const CosmicRayShield& cosmicRayShield ) const{
//     _straws.clear();
//     for ( size_t i=0; i<_indices.size(); ++i ){
//       StrawIndex idx = _indices[i];
//       const Straw* straw =  &cosmicRayShield.getStraw(idx);
//       _straws.push_back(straw);
//       straw->fillPointers(cosmicRayShield);
//     }
//   }

} // namespace mu2e 
#endif
  
