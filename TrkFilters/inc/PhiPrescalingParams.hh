#ifndef TrkFilters_PhiPrescalingParams_hh
#define TrkFilters_PhiPrescalingParams_hh

struct PhiPrescalingParams {
  float    _amplitude;
  float    _frequency;
  float    _phase;

  PhiPrescalingParams():
    _amplitude(0),
    _frequency(0),
    _phase(0) {}
  
  PhiPrescalingParams(float Amplitude, float Frequency, float Phase):
    _amplitude(Amplitude),
    _frequency(Frequency),
    _phase(Phase) {}
};


#endif /* TrkFilters_PhiPrescalingParams_hh */
