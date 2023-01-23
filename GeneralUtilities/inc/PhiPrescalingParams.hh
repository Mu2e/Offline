#ifndef GeneralUtilities_PhiPrescalingParams_hh
#define GeneralUtilities_PhiPrescalingParams_hh

#include "fhiclcpp/types/OptionalAtom.h"

namespace mu2e {
  struct PhiPrescalingParams {
    float    _amplitude;
    float    _frequency;
    float    _phase;

    struct Config {
      fhicl::OptionalAtom<float> amplitude{fhicl::Name("amplitude"), fhicl::Comment("amplitude")};
      fhicl::OptionalAtom<float> frequency{fhicl::Name("frequency"), fhicl::Comment("frequency")};
      fhicl::OptionalAtom<float> phase    {fhicl::Name("phase"),     fhicl::Comment("phase")};
    };

    PhiPrescalingParams():
      _amplitude(0),
      _frequency(0),
      _phase(0) {}

    PhiPrescalingParams(float Amplitude, float Frequency, float Phase):
      _amplitude(Amplitude),
      _frequency(Frequency),
      _phase    (Phase) {}

    PhiPrescalingParams(const Config& conf){
      float x;
      if (conf.amplitude(x)) _amplitude = x;
      if (conf.frequency(x)) _frequency = x;
      if (conf.phase(x))     _phase     = x;
    }
  };

}
#endif /* GeneralUtilities_PhiPrescalingParams_hh */
