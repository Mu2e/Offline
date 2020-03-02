// Simple approximations available for DIO spectrum.
//
// $Id: SimpleSpectrum.cc,v 1.6 2014/02/25 17:14:10 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/02/25 17:14:10 $
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
    if ( _approx > Spectrum::FlatTrunc && _phy.getCzarneckiCoefficients().empty() )
      throw cet::exception("EmptyCoefficients") <<
        " No Czarnecki coefficients available!  Cannot use this approximation.\n" ;
  }

  SimpleSpectrum::~SimpleSpectrum() {
  }

  double SimpleSpectrum::getWeight(double E) const {

    double weight(0.);
    if      ( _approx == Spectrum::Flat        ) weight = getFlat     (E, _phy);
    else if ( _approx == Spectrum::FlatTrunc   ) weight = getFlatTrunc(E, _phy);
    else if ( _approx == Spectrum::Pol5        ) weight = getPol5     (E, _phy);
    else if ( _approx == Spectrum::Pol58       ) weight = getPol58    (E, _phy);

    return weight;      

  }

  //========================================
  // Simple approximations below
  //========================================

  double SimpleSpectrum::getFlat(const double e, const PhysicsParams& phy ) {
    return  1.; 
  }

  double SimpleSpectrum::getFlatTrunc(const double e, const PhysicsParams& phy ) {
    return ( e > phy.getEndpointEnergy() ) ? 0. : 1.; 
  }

  double SimpleSpectrum::getPol5(const double e, const PhysicsParams& phy ) {

    const double delta = phy.getMuonEnergy() - e - cet::pow<2>( e )/(2*phy.getAtomicMass());
    
    return (e > phy.getEndpointEnergy() ) ? 0. : phy.getCzarneckiCoefficient()*cet::pow<5>( delta );

  }

  double SimpleSpectrum::getPol58(const double e, const PhysicsParams& phy ) {

    const double delta = phy.getMuonEnergy() - e - cet::pow<2>( e )/(2*phy.getAtomicMass());
      
    if ( e > phy.getEndpointEnergy() ) return 0.;
    
    const auto & coeffs = phy.getCzarneckiCoefficients();

    double prob(0.);
    double power = cet::pow<5>( delta );
    for ( size_t i=0; i < coeffs.size() ; i++ ) {
      if( i > 0 ) power *= delta; 
      prob  += coeffs.at(i)*power;
    }

    return prob;

  }

}
