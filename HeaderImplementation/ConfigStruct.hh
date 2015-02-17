#ifndef ConfigStruct_hh
#define ConfigStruct_hh

#include "TMath.h"

// This struct contains all parameters which remain constant throughout the simulation
struct ConfigStruct{
    const Double_t _shapingTime; // Shaping time (in units of ns)
    const Int_t _numSamplesPerHit; // Number of samples measured per hit
    const Double_t _adcError; // Assumes constant error for all adc measurements (in units of bits)
    const Double_t _measurementFrequency; // Sample frequency of adc values (in units of ns)
    const Double_t _truncationLevel; // Level of truncation of waveform (in units of bits)
    const Double_t _defaultPedestal; // Count value corresponding to the default pedestal (in units of bits)


    // Precomputed constants
    const Double_t _bits2scalingFactor; // approximately 67.96
    const Double_t _scalingFactor2bits; // approximately 0.0147
    const Double_t _hitPeriod;



    ConfigStruct() : _shapingTime(25.0), 
                     _numSamplesPerHit(8), 
                     _measurementFrequency(20.0), 
                     _adcError(3.0), 
                     _truncationLevel(1023.0),
                     _defaultPedestal(64.0),
                     /**_bits2scalingFactor(_shapingTime * TMath::E()),
                     _scalingFactor2bits(1.0 / _bits2scalingFactor),
                     _hitPeriod((_numSamplesPerHit - 1.0) * _measurementFrequency)**/
                     _bits2scalingFactor(25.0 * TMath::E()),
                     _scalingFactor2bits(1.0 / 25.0 / TMath::E()),
                     _hitPeriod(14.0)
                     {} 




}; 
#endif
