#ifndef Mu2eUtilities_MueXSpectrum_hh
#define Mu2eUtilities_MueXSpectrum_hh

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>

namespace mu2e {

  double f(double E, void *p);

  class MueXSpectrum {

  public:
    
    struct Params_t {
      double eMax;
      double me;
    } _par;
    
    MueXSpectrum(double maxEnergy, double bin);
    
    ~MueXSpectrum(){}
 
    double getWeight                     (double E) const;
    double getCorrectedMueXSpectrum(double e) const ;
    double evalIntegral                  (double de);
    static double  f                  (double E, void *p);

    void   setSpectrum   (int SpectrumType) { _spectrumType = SpectrumType; }
  
  private:

    double             _bin;
    int                _spectrumType;   
    double             _eMax;           
    double             _me;		// electron mass

    int                _nbins;
    double             _integral;      // over n-1 bins...
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_MueXSpectrum_hh */

  
