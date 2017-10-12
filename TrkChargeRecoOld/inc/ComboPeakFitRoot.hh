
#ifndef TrkChargeRecoOld_ComboPeakFitRoot_hh
#define TrkChargeRecoOld_ComboPeakFitRoot_hh
// base class for controlling peak fits
//
#include "TrkChargeRecoOld/inc/PeakFitRoot.hh"
#include "TGraphErrors.h"
#include <vector>

// base class for controlling the peak fit.
namespace mu2e {

  namespace TrkChargeRecoOld {

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
		ComboPeakFitRoot(StrawElectronics const& strawele, FitConfig const& config,FitType const& fittype, std::string fitoptions="QNEX0S");

		struct peakResult
		{
			Float_t _peakTime;  //time of peak relative to 140.0 ns interval
			Float_t _peakHeight; //in units of bits

			//Default constructor  
			peakResult() : _peakTime(0.0), _peakHeight(0.0){};

			//True constructor
			peakResult(Float_t peakTime, Float_t peakHeight) : _peakTime(peakTime), _peakHeight(peakHeight){};
    	};

    	//This is object will be filled by the findPeaks and addEarlyPeak methods 
    	typedef std::vector<peakResult> peakResultVector;

		bool hasEarlyCharge(const peakResultVector &initialGuess) const;

		bool hasLateCharge(const peakResultVector &initialGuess) const;
	
		void addEarlyPeak(const TGraphErrors &gr, peakResultVector &initialGuess) const;

		void findPeaks(const TGraphErrors &gr, peakResultVector &initialGuess, const double sigma = 3.0) const;
    };
  }
}
#endif
