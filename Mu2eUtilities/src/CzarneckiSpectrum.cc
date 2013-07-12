//
// Read the table with data about DIO spectrum, and
// merge the spectrum with the analytic expression
// in the endpoint region taken from Czarnecki spectrum
// Czarneckki et al 10.1103/PhysRevD.84.013006
//
// $Id: CzarneckiSpectrum.cc,v 1.8 2013/07/12 17:17:38 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/12 17:17:38 $
//

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"

// Framework includes
#include "cetlib/pow.h"

// C++ includes
#include <fstream>
#include <iostream>
#include <vector>

namespace mu2e {

  CzarneckiSpectrum::CzarneckiSpectrum()
  {
    readTable();
    checkTable();
  }

  double CzarneckiSpectrum::getWeight(double E) {

    std::vector<SpectrumValue>::iterator it = _table.begin();
    while ((it != _table.end()) && (it->first > E)) it++;

    if (it == _table.end()) return 0;
    
    double weight (0.);

    if (it->first <= E) { 
      if ( it==_table.begin() || (it+1)==_table.end()) {
	weight = it->second;
      } else {
        //        std::cout << "Interpolating" << std::endl;
        weight = interpolate(E, (it+1)->first, (it+1)->second,
                             it->first, it->second,
                             (it-1)->first, (it-1)->second);

        if ( weight < 0 ) weight = interpolateE5 ( E, *it );
        
      }
    }

    return weight;
  }

  void CzarneckiSpectrum::readTable() {

    ConfigFileLookupPolicy findConfig;

    GlobalConstantsHandle<PhysicsParams> phy;

    std::string filename = findConfig("ConditionsService/data/czarnecki_"+phy->getStoppingTarget()+".tbl");

    std::fstream intable(filename.c_str(),std::ios::in);
    if (!(intable.is_open())) {
      throw cet::exception("ProductNotFound")
        << "No Tabulated spectrum table file found";
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (!(intable.eof())) {
	_table.emplace_back( en, prob );
      }
    }

    // If the highest-energy entry does not have a weight of 0
    // insert one entry with 0, spaced equidistantly wrt to the highest-energy pt.
    if ( _table.begin()->second != 0. ) {
      auto it0 = _table.begin();
      auto it1 = it0+1;
      double spacing = it0->first - it1->first;
      _table.emplace( it0, it0->first + spacing, 0. );
    }

  }
  
  /*  double CzarneckiSpectrum::FitCzarnecki(double E) {
      
  double delta = 105.194 - E - E*E/(2*25133);
  
  double valueIs = (8.6434e-17)*pow(delta,5) + (1.16874e-17)*pow(delta,6) 
  - (1.87828e-19)*pow(delta,7) + (9.16327e-20)*pow(delta,8);
  
  return valueIs;
  
  }*/ 
  //Maybe in a later step we might want to use the fit function (valid from 85 MeV on)

  double CzarneckiSpectrum::interpolate(double E, double e1, double p1,
                                        double e2, double p2, double e3, double p3) {
    
    double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }
  
  double CzarneckiSpectrum::interpolateE5(double E, const SpectrumValue& val ) {
    
    GlobalConstantsHandle<PhysicsParams> phy;

    double b = val.second/cet::pow<5>( phy->getMuonEnergy()-val.first-cet::square(val.first)/(2*phy->getAtomicMass()));
    if ( phy->getEndpointEnergy() - E < 0 ) return 0.;
    else return b*cet::pow<5>( phy->getMuonEnergy() - E - cet::square(E)/(2*phy->getAtomicMass()) );
  }

}

