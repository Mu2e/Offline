#ifndef Mu2eUtilities_ConversionSpectrum_hh
#define Mu2eUtilities_ConversionSpectrum_hh

// Mu2e includes
// #include "Mu2eUtilities/inc/Table.hh"

// CLHEP includes
#include "CLHEP/Vector/LorentzVector.h"

// C++ includes
#include <utility>


/*namespace CLHEP {
  class RandFlat;
  }*/
namespace mu2e {


  double my_f(double E, void *p);


  class ConversionSpectrum {

  public:
    
    struct Params_t {
      double eMax;
      double alpha;
      double me;
    } _par;
    


    ConversionSpectrum(double maxEnergy, double bin, int RadCorrected = 0);
    
    ~ConversionSpectrum(){}
 
    double getWeight                     (double E) const;
    double getCorrectedConversionSpectrum(double e) const ;
    double evalIntegral                  (double de);
    static double  my_f                  (double E, void *p);

    void   setSpectrum   (int SpectrumType) { _spectrumType = SpectrumType; }
  
  private:

    double             _bin;
    int                _spectrumType;   // 0:delta function ; 1: rad corrected
    double             _eMax;           // max electron/positron energy
    double             _me;		// electron mass

    int                _nbins;
    double             _integral;      // over n-1 bins...
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_ConversionSpectrum_hh */

  
