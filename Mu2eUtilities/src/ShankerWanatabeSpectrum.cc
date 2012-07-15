//
// Read the table with wanatabe data about DIO spectrum, and
// merge the spectrum with the corrected Shanker analytic expression
// after the data endpoint.
//
// $Id: ShankerWanatabeSpectrum.cc,v 1.16 2012/07/15 22:06:18 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:18 $
//

#include "Mu2eUtilities/inc/ShankerWanatabeSpectrum.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "cetlib/pow.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

using cet::cube;
using cet::pow;
using cet::square;

namespace mu2e {

  ShankerWanatabeSpectrum::ShankerWanatabeSpectrum(int atomicZ, double mumass, double emass):
  //atomic number of the foil material
    _znum ( atomicZ ),
    _mumass ( mumass ),
    _emass ( emass )
  {
    readWanatabeTable();
    checkTable();
    //cout << "interpulation test f(x) = x^2+3x-5. f(3) = " << Interpulate(3,0,-5,2,5,10,125) << endl;
      _norm = _wanaEndPointVal / evaluateShanker(_wanaEndPoint);
  }

  ShankerWanatabeSpectrum::~ShankerWanatabeSpectrum()
  {
  }

  double ShankerWanatabeSpectrum::operator()(double E) {

    double value=0;
    if (E<_wanaEndPoint) {
      value = evaluateWanatabe(E);
    }
    if (E>=_wanaEndPoint) {
      value = _norm * evaluateShanker(E);
      if (E == _wanaEndPoint) cout << "Value at merging point is " << value << endl;
    }
    return value;

  }

  void ShankerWanatabeSpectrum::readWanatabeTable() {

    ConfigFileLookupPolicy findConfig;

    fstream intable(findConfig("ConditionsService/data/wanatabe.tbl").c_str(),ios::in);
    if (!(intable.is_open())) {
      throw cet::exception("ProductNotFound")
        << "No Wanatabe spectrum table file found";
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (en!=0&&prob!=0) {
	Value valueToAdd;
	valueToAdd.energy = en;
	valueToAdd.weight = prob;
        _wanatable.push_back(valueToAdd);
      }
    }

    _wanaEndPoint = _wanatable.front().energy;
    _wanaEndPointVal = _wanatable.front().weight;

  }

  void ShankerWanatabeSpectrum::checkTable() {

    double valueToCompare = (_wanatable.at(0).energy) + 1e9; 
    //order check
    for ( vector<Value>::iterator it =  _wanatable.begin(); it != _wanatable.end(); ++it) {
      if (it->energy >= valueToCompare) {
      throw cet::exception("Format")
        << "Wrong value in the wanatabe table: " << it->energy;
      }
    }
    //    unsigned tablesize = _table.size();

  }


  double ShankerWanatabeSpectrum::evaluateShanker(double E) {

    double BindEnergy = 13.6 * ( _mumass / _emass ) * _znum * _znum / 1e6;
    double atMassToMev = CLHEP::amu_c2;

    double AlAtWeight = 26.9815;
    double AtomicWeightMev = AlAtWeight * atMassToMev;

    double deltaPrimeMax = _mumass - BindEnergy;
    double muEndPoint = deltaPrimeMax - ((deltaPrimeMax*deltaPrimeMax)/(2*AtomicWeightMev));

    if (E > muEndPoint) return 0;

    double deltaOne = _mumass - BindEnergy - E - ((E*E)/(2*AtomicWeightMev));

    double ShD = 7.21861e-4*cube((double)_znum)-2.5289e-2*square((double)_znum)+0.388249*_znum-1.98475;
    double ShE = 1.76407e-3*cube((double)_znum)-5.19805e-2*square((double)_znum)+0.736126*_znum-3.69662;
    double ShF = 5.41126e-3*cube((double)_znum)-0.165584*square((double)_znum)+1.9329*_znum-6.7013;

    double shterm1 = E*E / (_mumass*_mumass);
    double shterm2 = pow<5>(deltaOne/_mumass);
    double shterm4 = ShE * deltaOne / _mumass;
    double shterm5 = ShF * (deltaPrimeMax - E);

    double TotShanker = shterm1 * shterm2 * (ShD+shterm4+shterm5);

    return TotShanker;

  }

  double ShankerWanatabeSpectrum::evaluateWanatabe(double E) {

    vector<Value>::iterator it = _wanatable.begin();
    while ((it != _wanatable.end()) &&(E < it->energy-0.0049)) {
      it++;
    }

    if (it == _wanatable.end()) {
      return 0;
    }

    if (it->energy <= E + 0.0049 || it->energy >= E - 0.0049 ) { //tollerance of 0.049 MeV
      return it->weight;
    } else {
      return interpulate(E, (it+1)->energy, (it+1)->weight,
                         it->energy, it->weight,
                         (it-1)->energy, (it-1)->weight);
    }

    return 0;

  }

  double ShankerWanatabeSpectrum::interpulate(double E, double e1, double p1,
                                              double e2, double p2, double e3, double p3) {

    double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }


}

