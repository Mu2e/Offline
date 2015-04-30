
#ifndef TrkChargeReco_ComboPeakFitRoot_hh
#define TrkChargeReco_ComboPeakFitRoot_hh
// base class for controlling peak fits
//
#include "TrkChargeReco/inc/PeakFitRoot.hh"
#include "TGraphErrors.h"
#include <vector>

// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeReco {

    //  PeakFitRoot
    //  Virtual class providing structure for extracting charge from ADC waveforms
    class ComboPeakFitRoot : public PeakFitRoot  {
      public:
	
	// extract peak information from adc waveform data.  1 waveform generates 1 peak fit.
	// The default implementation simply sums the ADC data after subtracting pedesdal
	virtual void process(StrawElectronics::ADCWaveform const& adcData, PeakFitParams & fit) const;

	// Destructor
	virtual ~ComboPeakFitRoot(){}

	// ComboPeakFitRoot normal constructor with ConfigStruct initilization parameters
	ComboPeakFitRoot(StrawElectronics const& strawele, FitConfig const& config,std::string fitoptions="QNEX0S");

	struct resultantPeakData
	{
		Float_t _peakTime;  //time of peak relative to 140.0 ns interval
		Float_t _peakHeight; //in units of bits

		//Default constructor  
		resultantPeakData() : _peakTime(0.0), _peakHeight(0.0){};

		//True constructor
		resultantPeakData(Float_t peakTime, Float_t peakHeight) : _peakTime(peakTime), _peakHeight(peakHeight){};
    };

    //This is object will be filled by the findPeaks and addEarlyPeak methods 
    	typedef std::vector<resultantPeakData> resultantHitData;

	bool hasEarlyCharge(const resultantHitData &initialGuess) const;

	bool hasLateCharge(const resultantHitData &initialGuess) const;

	bool hasFloatingPedestal(StrawElectronics::ADCWaveform const& adcData) const;	

	void addEarlyPeak(const TGraphErrors &gr, resultantHitData &initialGuess) const;

	void findPeaks(const TGraphErrors &gr, resultantHitData &initialGuess, const double sigma = 3.0) const;
    };
  }
}
#endif
