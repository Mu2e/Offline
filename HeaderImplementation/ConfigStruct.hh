#ifndef ConfigStruct_hh
#define ConfigStruct_hh

// This struct contains all parameters which remain constant throughout the simulation
struct ConfigStruct{
    const Double_t _shapingTime; // Shaping time (in units of ns)
    const Int_t _numSamplesPerHit; // Number of samples measured per hit
    const Double_t _adcError; // Assumes constant error for all adc measurements (in units of bits)
    const Double_t _measurementFrequency; // Sample frequency of adc values (in units of ns)
    const Double_t _truncationLevel; // Level of truncation of waveform (in units of bits)
    const Double_t _defaultPedestal; // Count value corresponding to the default pedestal (in units of bits)

    ConfigStruct() : _shapingTime(25.0), 
                     _numSamplesPerHit(8), 
                     _measurementFrequency(20.0), 
                     _adcError(3.0), 
                     _truncationLevel(1023.0),
                     _defaultPedestal(64.0){} 
}; 
#endif
