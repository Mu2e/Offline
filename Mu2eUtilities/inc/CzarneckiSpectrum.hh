#ifndef Mu2eUtilities_CzarneckiSpectrum_hh
#define Mu2eUtilities_CzarneckiSpectrum_hh
//
// Read Czarnecki DIO spectrum from a table and merge it
// with the spectrum coming from the endopoint region formula

// $Id: CzarneckiSpectrum.hh,v 1.4 2012/05/08 19:09:59 onoratog Exp $
// $Author: onoratog $
// $Date: 2012/05/08 19:09:59 $
//
// Original Author: Gianni Onorato
//

// C++ includes
#include <utility>
#include <vector>

// Framework includes
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes
#include "Mu2eUtilities/inc/DIOBase.hh"

//CLHEP includes
#include "CLHEP/Random/RandGeneral.h"

namespace mu2e {

  class CzarneckiSpectrum {

    struct Value{
      double energy;
      double weight;
    };
      
      
  public:
    
    CzarneckiSpectrum(int atomicZ);
    
    ~CzarneckiSpectrum();

    double operator()(double E);

  private:

    int _znum;

    std::vector<Value> _table;

    void readTable();

    void checkTable();

    //    double FitCzarnecki(double E); maybe we'll use it later

    double interpulate(double E, double e1, double p1,
                       double e2, double p2, double e3, double p3);  

  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_CzarneckiSpectrum_hh */

