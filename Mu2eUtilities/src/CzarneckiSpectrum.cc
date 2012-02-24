//
// Read the table with data about DIO spectrum, and
// merge the spectrum with the analytic expression
// in the endpoint region taken from Czarnecki.
//
// $Id: CzarneckiSpectrum.cc,v 1.3 2012/02/24 20:05:52 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/02/24 20:05:52 $
//

#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "Mu2eUtilities/inc/ConfigFileLookupPolicy.hh"
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

  CzarneckiSpectrum::CzarneckiSpectrum(int atomicZ):
  //atomic number of the foil material
    _znum ( atomicZ )
  {
    readTable();
  }

  CzarneckiSpectrum::~CzarneckiSpectrum()
  {
  }

  double CzarneckiSpectrum::operator[](double E) {

    vector<pair<double, double> >::iterator it = _table.begin();
    //    cout << "Searching for " << E << endl;
    while ((E < it->first-0.0049) && it != _table.end()) {
      //    cout << "In the table I have " << it->first << endl;
      it++;
    }

    if (it == _table.end()) {
      return 0;
    } 

    if (it->first <= E + 0.0049 || it->first >= E - 0.0049 ) { //tollerance of 0.0049 MeV
      //  cout << "And so I assign " << it->second << endl;
      return it->second;
    } else {
      // cout << "Interpulating" << endl;
      return interpulate(E, (it+1)->first, (it+1)->second,
                         it->first, it->second,
                         (it-1)->first, (it-1)->second);
    }

    //    if (it->first < E) {
    //  cout << "Assignin through fit" << endl;
    //  return FitCzarnecki(E);
    // }

    return 0;

  }

  void CzarneckiSpectrum::readTable() {

    ConfigFileLookupPolicy findConfig;

    fstream intable(findConfig("ConditionsService/data/czarnecki.tbl").c_str(),ios::in);
    if (!(intable.is_open())) {
      throw cet::exception("ProductNotFound")
        << "No Tabulated spectrum table file found";
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (en!=0&&prob!=0)
        _table.push_back(pair<double,double>(en,prob));
    }

  }

  /*  double CzarneckiSpectrum::FitCzarnecki(double E) {
      
  double delta = 105.194 - E - E*E/(2*25133);
  
  double valueIs = (8.6434e-17)*pow(delta,5) + (1.16874e-17)*pow(delta,6) 
  - (1.87828e-19)*pow(delta,7) + (9.16327e-20)*pow(delta,8);
  
  return valueIs;
  
  }*/ 
  //Maybe in a later step we might want to use the fit function (valid from 85 MeV on)

  double CzarneckiSpectrum::interpulate(double E, double e1, double p1,
                                        double e2, double p2, double e3, double p3) {
    
    double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }
  
}

