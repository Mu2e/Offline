#ifndef DiagPeakFit_hh 
#define DiagPeakFit_hh

#include "TrkChargeReco/inc/ComboPeakFit.hh"

// Contains PeakFit classes which should be used for diagnostic purposes

namespace mu2e {

namespace TrkChargeReco {

class SinglePeakWithoutTruncFit : public SinglePeakFit{
	public:
		// SinglePeakFit normal constructor with ConfigStruct initilization parameters
		SinglePeakWithoutTruncFit(const ConfigStruct &initParams);
};

class SinglePeakFloatingPedestalWithoutTruncFit : public SinglePeakFloatingPedestalFit{
	public:
		// FindSinglePeak normal constructor with ConfigStruct initilization parameters
		SinglePeakFloatingPedestalWithoutTruncFit(const ConfigStruct &initParams);
};

class  LXPeakVsSinglePeakFit : public PeakFitRootBase{
	public:
		LXPeakVsSinglePeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams){};

		virtual void process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result);
};
}
}
#endif
