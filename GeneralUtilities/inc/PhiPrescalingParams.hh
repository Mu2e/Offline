#ifndef GeneralUtilities_PhiPrescalingParams_hh
#define GeneralUtilities_PhiPrescalingParams_hh

#include "fhiclcpp/types/Atom.h"

namespace mu2e {
  struct PhiPrescalingParams {
    float    _amplitude;
    float    _frequency;
    float    _phase;

    struct Config {
      fhicl::Atom<float> amplitude{fhicl::Name("amplitude"), fhicl::Comment("amplitude")};
      fhicl::Atom<float> frequency{fhicl::Name("frequency"), fhicl::Comment("frequency")};
      fhicl::Atom<float> phase    {fhicl::Name("phase"),     fhicl::Comment("phase")};
    };
    
    PhiPrescalingParams():
      _amplitude(0),
      _frequency(0),
      _phase(0) {}
    
    PhiPrescalingParams(float Amplitude, float Frequency, float Phase):
      _amplitude(Amplitude),
      _frequency(Frequency),
      _phase    (Phase) {}
    
    PhiPrescalingParams(const Config& conf): 
      _amplitude(conf.amplitude()),
      _frequency(conf.frequency()),
      _phase    (conf.phase()) {}
  };

}
#endif /* GeneralUtilities_PhiPrescalingParams_hh */
