#ifndef FindPeakBase_hh
#define FindPeakBase_hh

#include "TMath.h"
#include "ConfigStruct.hh"

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

/** \class FindPeakBase
 * Virtual class providing structure for FindSinglePeak, FindDoublePeak, FindMutiplePeaks, etc. 
 */
 class FindPeakBase{
  public:
    
    // Fills result using adc waveform data
    // NOTE : This function may begin with peak data provided in result which is replaced
    virtual void process(const adcWaveform adcData, resultantHitData &result) = 0;

    // Destructor
    virtual ~FindPeakBase(){}

    // FindPeakBase normal constructor with ConfigStruct initilization parameters
    FindPeakBase(const ConfigStruct &initParams) : _initParams(initParams){}
  protected:

    ConfigStruct _initParams; 
};
#endif