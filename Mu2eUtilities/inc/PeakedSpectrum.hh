#ifndef Mu2eUtilities_PeakedSpectrum_hh
#define Mu2eUtilities_PeakedSpectrum_hh

// C++ includes
#include <utility>

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"


namespace CLHEP { class RandFlat;     }
namespace fhicl { class ParameterSet; }

namespace mu2e {

  class PeakedSpectrum {

  public:

    enum enum_type    { kPowerLaw, kTriangleFlat };

    // random number generators ar owned by the callers, no memory cleanup needed
    PeakedSpectrum(const fhicl::ParameterSet& psphys, CLHEP::RandFlat& randFlat);

    ~PeakedSpectrum(){}

    double getWeight   (double x) const;
    double getMax      () const;

    void   setSpectrum   (enum_type    spectrum  ) { _spectrum   = spectrum;   }
    double fire() const;


  private:
    enum { kMaxParameters = 10 };

    CLHEP::RandFlat&   _rnFlat;
    int                _verbose;

    enum_type          _spectrum;
    double             _norm;
    double             _xmax;
    double             _xmin;
    double             _paramsD[kMaxParameters];
    int                _paramsI[kMaxParameters];
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_PeakedSpectrum_hh */
