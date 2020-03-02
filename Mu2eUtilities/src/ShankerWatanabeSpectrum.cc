//
// Read the table with Watanabe data about DIO spectrum, and
// merge the spectrum with the corrected Shanker analytic expression
// after the data endpoint.
//
// $Id: ShankerWatanabeSpectrum.cc,v 1.6 2014/05/01 18:12:26 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2014/05/01 18:12:26 $
//

// Mu2e includes
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/PhysicsParams.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/ShankerWatanabeSpectrum.hh"

// CLHEP includes
#include "CLHEP/Units/PhysicalConstants.h"

// Framework includes
#include "cetlib/pow.h"

// C++ includes
#include <fstream>
#include <iostream>

using namespace std;

namespace mu2e {

  ShankerWatanabeSpectrum::ShankerWatanabeSpectrum() :
    _table ( loadTable<2,false>( ConfigFileLookupPolicy()("ConditionsService/data/watanabe.tbl" ) ) )
  {
    _wanaEndPoint    = _table(0,0); 
    _wanaEndPointVal = _table(0,1);
    _norm     = _wanaEndPointVal / evaluateShanker( _wanaEndPoint );
  }


  double ShankerWatanabeSpectrum::getWeight(double E) const {

    double value(0.);
    if ( E < _wanaEndPoint ) {
      value = evaluateWatanabe(E);
    }
    if (E>=_wanaEndPoint) {
      value = _norm * evaluateShanker(E);
      if (E == _wanaEndPoint) cout << "Value at merging point is " << value << endl;
    }
    return value;

  }

  double ShankerWatanabeSpectrum::evaluateShanker(double E) const {

    GlobalConstantsHandle<PhysicsParams>     phy;
    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Shanker uses approx. binding energy
    const double bindEnergy = phy->getApproxEb();

    // pick up particle mass
    const double mumass = pdt->particle(PDGCode::mu_minus).ref().mass().value();

    const double deltaPrimeMax = mumass - bindEnergy;
    const double muEndPoint    = deltaPrimeMax - cet::square(deltaPrimeMax)/(2*phy->getAtomicMass());
    
    if (E > muEndPoint) return 0;
    
    const double delta1 = mumass - bindEnergy - E - cet::square(E)/(2*phy->getAtomicMass());

    double shD(0.), shE(0.), shF(0.);
    unsigned zpower (1);
    
    for ( size_t i(0); i < phy->getShankerNcoeffs() ; i++ ) {
      shD    += phy->getShankerDcoefficients().at(i)*zpower;
      shE    += phy->getShankerEcoefficients().at(i)*zpower;
      shF    += phy->getShankerFcoefficients().at(i)*zpower;
      zpower *= phy->getAtomicNumber();
    }

    const double shterm1 = cet::square(E/mumass);
    const double shterm2 = cet::pow<5>(delta1/mumass);
    const double shterm4 = shE * delta1 / mumass;
    const double shterm5 = shF * (deltaPrimeMax - E);

    return shterm1 * shterm2 * (shD+shterm4+shterm5);

  }

  double ShankerWatanabeSpectrum::evaluateWatanabe(double E) const {

    double weight (0.);

    const unsigned iRow = _table.getLowerBoundRowIndex( E );
    if ( iRow == 0 || iRow == _table.getNrows() ) return weight;

    auto const & row        = _table.getRow( iRow   );
    auto const & row_before = _table.getRow( iRow-1 );
    auto const & row_after  = _table.getRow( iRow+1 );

    return interpolate(E, row_before, row, row_after );

  }


  double ShankerWatanabeSpectrum::interpolate(const double E, 
                                              const TableRow<2>& row_after,
                                              const TableRow<2>& row,
                                              const TableRow<2>& row_before ) const {
    
    const double e1(  row_after.first );  const double p1(  row_after.second.at(0) );
    const double e2(        row.first );  const double p2(        row.second.at(0) );
    const double e3( row_before.first );  const double p3( row_before.second.at(0) );
    
    const double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;
    
    const double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    const double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    const double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                      e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }

}

