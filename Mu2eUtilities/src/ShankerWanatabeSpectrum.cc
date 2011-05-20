//
// Read the table with wanatabe data about DIO spectrum, and
// merge the spectrum with the corrected Shanker analytic expression
// after the data endpoint.
//
// $Id: ShankerWanatabeSpectrum.cc,v 1.9 2011/05/20 22:39:28 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 22:39:28 $
//

#include "Mu2eUtilities/inc/ShankerWanatabeSpectrum.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "cetlib/pow.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

using cet::cube;
using cet::pow;
using cet::square;

namespace mu2e {

  ShankerWanatabeSpectrum::ShankerWanatabeSpectrum(int atomicZ):
  //atomic number of the foil material
    _Znum ( atomicZ )
  {
    ReadWanatabeTable();
    //cout << "interpulation test f(x) = x^2+3x-5. f(3) = " << Interpulate(3,0,-5,2,5,10,125) << endl;
      _norm = _WanaEndPointVal / EvaluateShanker(_WanaEndPoint);
  }

  ShankerWanatabeSpectrum::~ShankerWanatabeSpectrum()
  {
  }

  double ShankerWanatabeSpectrum::operator[](double E) {

    double value=0;
    if (E<_WanaEndPoint) {
      value = EvaluateWanatabe(E);
    }
    if (E>=_WanaEndPoint) {
      value = _norm * EvaluateShanker(E);
      if (E == _WanaEndPoint) cout << "Value at merging point is " << value << endl;
    }
    return value;

  }

  void ShankerWanatabeSpectrum::ReadWanatabeTable() {

    ConfigFileLookupPolicy findConfig;

    fstream intable(findConfig("ConditionsService/data/wanatabe.tbl").c_str(),ios::in);
    if (!(intable.is_open())) {
      throw cet::exception("ProductNotFound")
        << "No Wanatabe spectrum table file found";
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (en!=0&&prob!=0)
        _wanatable.push_back(pair<double,double>(en,prob));
    }

    _WanaEndPoint = _wanatable.front().first;
    _WanaEndPointVal = _wanatable.front().second;

  }

  double ShankerWanatabeSpectrum::EvaluateShanker(double E) {

    ConditionsHandle<ParticleDataTable> pdt("ignored");
    const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus).ref();
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus).ref();
    double mumass = mu_data.mass().value();
    double emass = e_data.mass().value();
    double BindEnergy = 13.6 * ( mumass / emass ) * _Znum * _Znum / 1e6;
    double atMassToMev = CLHEP::amu_c2;

    double AlAtWeight = 26.9815;
    double AtomicWeightMev = AlAtWeight * atMassToMev;

    double deltaPrimeMax = mumass - BindEnergy;
    double muEndPoint = deltaPrimeMax - ((deltaPrimeMax*deltaPrimeMax)/(2*AtomicWeightMev));

    if (E > muEndPoint) return 0;

    double deltaOne = mumass - BindEnergy - E - ((E*E)/(2*AtomicWeightMev));

    double ShD = 7.21861e-4*cube((double)_Znum)-2.5289e-2*square((double)_Znum)+0.388249*_Znum-1.98475;
    double ShE = 1.76407e-3*cube((double)_Znum)-5.19805e-2*square((double)_Znum)+0.736126*_Znum-3.69662;
    double ShF = 5.41126e-3*cube((double)_Znum)-0.165584*square((double)_Znum)+1.9329*_Znum-6.7013;

    double shterm1 = E*E / (mumass*mumass);
    double shterm2 = pow<5>(deltaOne/mumass);
    double shterm4 = ShE * deltaOne / mumass;
    double shterm5 = ShF * (deltaPrimeMax - E);

    double TotShanker = shterm1 * shterm2 * (ShD+shterm4+shterm5);

    return TotShanker;

  }

  double ShankerWanatabeSpectrum::EvaluateWanatabe(double E) {

    vector<pair<double, double> >::iterator it = _wanatable.begin();
    while (E < it->first && it != _wanatable.end()) {
      it++;
    }

    if (it == _wanatable.end()) {
      return 0;
    }

    if (it->first == E) {
      return it->second;
    }

    if (it->first < E) {
      return Interpulate(E, (it+1)->first, (it+1)->second,
                         it->first, it->second,
                         (it-1)->first, (it-1)->second);
    }

    return 0;

  }

  double ShankerWanatabeSpectrum::Interpulate(double E, double e1, double p1,
                                              double e2, double p2, double e3, double p3) {

    double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }


}

