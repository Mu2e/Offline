// Simple approximations available for DIO spectrum.
//
// $Id: SimpleSpectrum.cc,v 1.3 2013/07/24 18:48:24 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/24 18:48:24 $
//

// Mu2e includes
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"

// Framework includes
#include "cetlib/pow.h"

using namespace std;

using cet::pow;
using cet::square;

namespace mu2e {

  SimpleSpectrum::SimpleSpectrum( Spectrum::enum_type approx, 
                                  const PhysicsParams& physicsParams ) :
    _approx ( approx ) ,  
    _phy    ( physicsParams )
  {
    if ( _approx > Spectrum::Flat && _phy.getCzarneckiCoefficients().empty() )
      throw cet::exception("EmptyCoefficients") <<
        " No Czarnecki coefficients available!  Cannot use this approximation.\n" ;
  }

  SimpleSpectrum::~SimpleSpectrum() {
  }

  double SimpleSpectrum::getWeight(double E) {

    double weight(0.);
    if      ( _approx == Spectrum::Flat  ) weight = getFlat (E);
    else if ( _approx == Spectrum::Pol5  ) weight = getPol5 (E);
    else if ( _approx == Spectrum::Pol58 ) weight = getPol58(E);

    return weight;      

  }

  //========================================
  // Simple approximations below
  //========================================

  double SimpleSpectrum::getFlat(const double e, const PhysicsParams& phy ) {
    return ( e > phy.getEndpointEnergy() ) ? 0. : 1.; 
  }

  double SimpleSpectrum::getPol5(const double e, const PhysicsParams& phy ) {

    const double delta = phy.getMuonEnergy() - e - cet::pow<2>( e )/(2*phy.getAtomicMass());
    
    if ( phy.getEndpointEnergy()-e < 0 ) return 0.;
    
    return phy.getCzarneckiCoefficient()*cet::pow<5>( delta );

  }

  double SimpleSpectrum::getPol58(const double e, const PhysicsParams& phy ) {

    const double delta = phy.getMuonEnergy() - e - cet::pow<2>( e )/(2*phy.getAtomicMass());
      
    if ( phy.getEndpointEnergy()-e < 0 ) return 0.;
    
    const auto & coeffs = phy.getCzarneckiCoefficients();

    double prob(0.);
    double power = cet::pow<5>( delta );
    for ( size_t i=0; i < coeffs.size() ; i++ ) {
      if( i > 0 ) power *= delta; 
      prob += coeffs.at(i)*power;
    }
    return prob;

  }

}
