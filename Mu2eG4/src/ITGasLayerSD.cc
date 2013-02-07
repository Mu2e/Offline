//
//
//  $Id: ITGasLayerSD.cc,v 1.15 2013/02/07 17:56:03 genser Exp $
//  $Author: genser $
//  $Date: 2013/02/07 17:56:03 $
//
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "Mu2eG4/inc/ITGasLayerSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"

#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  G4ThreeVector ITGasLayerSD::_mu2eDetCenter;

  ITGasLayerSD::ITGasLayerSD(G4String name, const SimpleConfig& config) :
                  Mu2eSensitiveDetector(name,config),
                  _superlayer(0),
                  _ring(0),
                  _nwires(0),
                  _Dphi(0)
  {
          art::ServiceHandle<GeometryService> geom;

          if ( !geom->hasElement<ITracker>() ) {
                  throw cet::exception("GEOM")
                  << "Expected I Trackers but found neither.\n";
          }
  }


  ITGasLayerSD::~ITGasLayerSD(){ }

} //namespace mu2e
