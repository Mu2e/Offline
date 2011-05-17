#ifndef SHANKERWANATABESPECTRUM_HH
#define SHANKERWANATABESPECTRUM_HH
//
// Read Wanatabe data about DIO spectrum from a table and merge it 
// with the spectrum coming from the Shanker formula

// $Id: ShankerWanatabeSpectrum.hh,v 1.2 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// 

// C++ incldues
#include <vector>
#include <utility>

// Framework includes
#include "art/Framework/Core/RandomNumberGeneratorService.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class ShankerWanatabeSpectrum {
    
  public:

    ShankerWanatabeSpectrum(int atomicZ);
    
    ~ShankerWanatabeSpectrum();
    
    double operator[](double E);

  private:

    int _Znum;

    double _WanaEndPoint, _WanaEndPointVal, _norm;

    std::vector<std::pair<double, double> > _wanatable;

    void ReadWanatabeTable();

    double EvaluateShanker(double E);

    double EvaluateWanatabe(double E);

    double Interpulate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);


  };

} // end of namespace mu2e

#endif

