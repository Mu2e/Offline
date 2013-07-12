//
// Read the table with Watanabe data about DIO spectrum, and
// merge the spectrum with the corrected Shanker analytic expression
// after the data endpoint.
//
// $Id: ShankerWatanabeSpectrum.cc,v 1.1 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//

// Mu2e includes
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MCDataProducts/inc/PDGCode.hh"
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

  ShankerWatanabeSpectrum::ShankerWatanabeSpectrum()
  {
    readTable();
    checkTable();
    //cout << "interpulation test f(x) = x^2+3x-5. f(3) = " << Interpulate(3,0,-5,2,5,10,125) << endl;
      _norm = _wanaEndPointVal / evaluateShanker(_wanaEndPoint);
  }


  double ShankerWatanabeSpectrum::getWeight(double E) {

    double value=0;
    if (E<_wanaEndPoint) {
      value = evaluateWatanabe(E);
    }
    if (E>=_wanaEndPoint) {
      value = _norm * evaluateShanker(E);
      if (E == _wanaEndPoint) cout << "Value at merging point is " << value << endl;
    }
    return value;

  }

  void ShankerWatanabeSpectrum::readTable() {

    ConfigFileLookupPolicy findConfig;

    fstream intable(findConfig("ConditionsService/data/wanatabe.tbl").c_str(),ios::in);
    if (!(intable.is_open())) {
      throw cet::exception("ProductNotFound")
        << "No Watanabe spectrum table file found";
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (en!=0&&prob!=0) {
        _table.emplace_back( en, prob );
      }
    }

    _wanaEndPoint    = _table.front().first;
    _wanaEndPointVal = _table.front().second;

  }

  double ShankerWatanabeSpectrum::evaluateShanker(double E) {

    GlobalConstantsHandle<PhysicsParams>     phy;
    GlobalConstantsHandle<ParticleDataTable> pdt;

    // Shanker uses approx. binding energy
    const double bindEnergy = phy->getApproxEb();

    // pick up particle mass
    const double mumass = pdt->particle(PDGCode::mu_minus).ref().mass().value();

    double deltaPrimeMax = mumass - bindEnergy;
    double muEndPoint = deltaPrimeMax - cet::square(deltaPrimeMax)/(2*phy->getAtomicMass());
    
    if (E > muEndPoint) return 0;
    
    double delta1 = mumass - bindEnergy - E - cet::square(E)/(2*phy->getAtomicMass());

    double shD(0.), shE(0.), shF(0.);
    unsigned zpower (1);
    
    for ( size_t i(0); i < phy->getShankerNcoeffs() ; i++ ) {
      shD    += phy->getShankerDcoefficients().at(i)*zpower;
      shE    += phy->getShankerEcoefficients().at(i)*zpower;
      shF    += phy->getShankerFcoefficients().at(i)*zpower;
      zpower *= phy->getAtomicNumber();
    }

    double shterm1 = cet::square(E/mumass);
    double shterm2 = cet::pow<5>(delta1/mumass);
    double shterm4 = shE * delta1 / mumass;
    double shterm5 = shF * (deltaPrimeMax - E);

    return shterm1 * shterm2 * (shD+shterm4+shterm5);

  }

  double ShankerWatanabeSpectrum::evaluateWatanabe(double E) {

    vector<SpectrumValue>::iterator it = _table.begin();
    while ((it != _table.end()) &&(E < it->first-_tolerance)) {
      it++;
    }

    if (it == _table.end()) {
      return 0;
    }

    if (it->first <= E + _tolerance || it->first >= E - _tolerance ) { //tolerance of 0.049 MeV
      return it->second;
    } else {
      return interpolate(E, (it+1)->first, (it+1)->second,
                         it->first, it->second,
                         (it-1)->first, (it-1)->second);
    }

    return 0;

  }

  double ShankerWatanabeSpectrum::interpolate(double E, double e1, double p1,
                                              double e2, double p2, double e3, double p3) {

    double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }


}

