//
// Generate an energy value for the DIO electrons, using Wanatabe
// data, merged to Shanker's formula near the endpoint
//
// $Id: DIOShankerWanatabe.cc,v 1.1 2011/02/28 16:18:50 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/02/28 16:18:50 $
//
// 

// C++ includes
#include <iostream>
#include <fstream>
#include <cmath>

// Framework includes
#include "FWCore/Utilities/interface/Exception.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOShankerWanatabe.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

using namespace std;

namespace mu2e {

  DIOShankerWanatabe::DIOShankerWanatabe(int atomicZ, double emin, double emax,
                                         edm::RandomNumberGeneratorService::base_engine_t& engine):
  //atomic number of the foil material
    _Znum ( atomicZ ),
  //limits on energy generation
    _emin ( emin ),
    _emax ( emax ),
    _randEnergy ( engine, &(ShWaSpectrum()[0]), _nBinsSpectrum )
  {
    if (_Znum!=13) {
      throw cms::Exception("GEOM")
        << "Foil material different from Alluminum.";
    }
  }

  DIOShankerWanatabe::~DIOShankerWanatabe()
  {
  }

  double DIOShankerWanatabe::fire() {

    double en = -1;
    do {
      en = _startpoint + (_endpoint-_startpoint)*_randEnergy.fire();
    } while (!(en >= _emin && en <= _emax));
    return en;
 
  }

  vector<double> DIOShankerWanatabe::ShWaSpectrum() {

    vector<double> spectrum;
    list<pair<double,double> > wanatable;
    list<pair<double,double> > shanktable;
    fstream intable("ConditionsService/data/wanatabe.tbl",ios::in);
    if (!(intable.is_open())) {
      cout << "No file found " << endl;
      //ADD EXCEPTION
      return spectrum;
    }
    double en, prob;
    while (!(intable.eof())) {
      intable >> en >> prob;
      if (en!=0&&prob!=0) 
        wanatable.push_back(pair<double,double>(en,prob));
    }
    AddShanker(wanatable.front(),shanktable);
    _nBinsSpectrum = wanatable.size()+shanktable.size()-1;
    list<pair<double,double> >::reverse_iterator wit = wanatable.rbegin();
    while (wit != wanatable.rend()) {
      double binVal = wit->second;
      spectrum.push_back(binVal);
      //     cout << "wata after build: " << wit->first << "  " << wit->second << endl;
      wit++;
    } 
    spectrum.pop_back(); 
    list<pair<double,double> >::reverse_iterator sit = shanktable.rbegin();
    while (sit != shanktable.rend()) {
      double binVal = sit->second;
      spectrum.push_back(binVal);
      //      cout << "shank after build: " << sit->first << "  " << sit->second << endl;
      sit++;
    } 
    _startpoint = wanatable.back().first;

    return spectrum;
  }
  

  void DIOShankerWanatabe::AddShanker(pair<double,double> wEnd, list<pair<double,double> > & shankList) {

    double NormVal = wEnd.second / EvalShanker(wEnd.first); 

    _endpoint = (double)(((int)(_endpoint*10+0.5))/10);

    //    cout << "endpoint: " << _endpoint << endl;

    double en = _endpoint;

    while (en>=wEnd.first) {
      double val = NormVal*EvalShanker(en);
      shankList.push_back(pair<double,double>(en,val));
      en = en - 0.1;
    }

  }

  double DIOShankerWanatabe::EvalShanker(double E) {

    ConditionsHandle<ParticleDataTable> pdt("ignored");
    const HepPDT::ParticleData& mu_data = pdt->particle(PDGCode::mu_minus);
    const HepPDT::ParticleData& e_data = pdt->particle(PDGCode::e_minus);
    double mumass = mu_data.mass().value();
    double emass = e_data.mass().value();
    double BindEnergy = 13.6 * ( mumass / emass ) * _Znum * _Znum / 1e6; 
    double atMassToMev = 931.0;

    double AlAtWeight = 26.9815;

    double AtomicWeightMev = AlAtWeight * atMassToMev;

    double deltaPrimeMax = mumass - BindEnergy;
    double muEndPoint = deltaPrimeMax - ((deltaPrimeMax*deltaPrimeMax)/(2*AtomicWeightMev));

    _endpoint = muEndPoint;

    double deltaOne = mumass - BindEnergy - E - ((E*E)/2*AtomicWeightMev);

    double ShD13 = 0.3745779;
    double ShE13 = 0.9639610;
    double ShF13 = 2.3311688;

    double shterm1 = E*E / mumass*mumass;
    double shterm2 = pow(deltaOne/mumass,5);
    double shterm4 = ShE13 * deltaOne / mumass;
    double shterm5 = ShF13 * (deltaPrimeMax - E);

    double TotShanker = shterm1 * shterm2 * (ShD13+shterm4+shterm5);

    return TotShanker;

  }

}

