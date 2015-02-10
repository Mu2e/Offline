#ifndef FindPeakBase_hh
#define FindPeakBase_hh

#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>
#include "config.hh"

typedef unsigned int * adcWaveform;

struct resultantPeakData
{
  Float_t _peakHeight; // in units of bits
  Float_t _peakTime; // time of peak relative to 140.0 ns interval

  // Default constructor - should probably be deleted
  resultantPeakData() : _peakTime(0.0), _peakHeight(0.0){};

  // True constructor
  resultantPeakData(Float_t peakTime, Float_t peakHeight) : _peakTime(peakTime), _peakHeight(peakHeight){};
};

// This is object top which will be filled by the process method 
typedef std::vector<resultantPeakData> resultantHitData;


// Virtual class providing structure for FindSinglePeak, FindDoublePeak, FindMutiplePeaks, etc. 
class FindPeakBase{
  public:
    
    // Fills result using adc waveform data
    // NOTE : This function may begin with peak data provided in result which is replaced
    virtual void process(const adcWaveform adcData, resultantHitData &result) = 0;

    // Destructor
    virtual ~FindPeakBase(){}

    // FindPeakBase normal constructor with configStruct initilization parameters
    FindPeakBase(const configStruct &initParams) : _initParams(initParams), 
                    _bits2scalingFactor(initParams._shapingTime * TMath::E()),
                    _scalingFactor2bits(1.0 / _bits2scalingFactor),
                    _hitPeriod(initParams._numSamplesPerHit * initParams._measurementFrequency - 1.0){}

  protected:

    configStruct _initParams; 

    // Precomputed constants
    const Double_t _bits2scalingFactor; // approximately 67.96
    const Double_t _scalingFactor2bits; // approximately 0.0147
    const Double_t _hitPeriod;
};
#endif