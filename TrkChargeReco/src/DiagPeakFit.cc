// MAKE SURE THAT initParams is getting passed everywhere correctly 

#include "TrkChargeReco/inc/FitModelRoot.hh"
#include "TrkChargeReco/inc/DiagPeakFit.hh"

// TODO : ADD CLASS LXPeakFitFloatingPedestal

namespace mu2e{

namespace TrkChargeReco{
// FindSinglePeakWithoutTrunc normal constructor with ConfigStruct initilization parameters
// NOTE THIS A HORRIBLE IDEA TO INHERIT THE CONSTRUCTOR OF SINGLE PEAK SINCE IT REPLACES _FITMODEL
SinglePeakWithoutTruncFit::SinglePeakWithoutTruncFit(const ConfigStruct &initParams) : SinglePeakFit(initParams)
{
  _fitModel = TF1("fitModel",FitModelRoot::singlePeak,0.0,_initParams._hitPeriod,3);
}

// FindSinglePeakFloatingPedestal normal constructor with ConfigStruct initilization parameters
SinglePeakFloatingPedestalWithoutTruncFit::SinglePeakFloatingPedestalWithoutTruncFit(const ConfigStruct &initParams) : SinglePeakFloatingPedestalFit(initParams)
{
  _fitModel = TF1("fitModel",FitModelRoot::singlePeakFloatingPedestal,0.0,_initParams._hitPeriod,4);
}

// Fills result using adc waveform data using either FindSinglePeak or LXPeakFit depending on the results of
// the peak search
void LXPeakVsSinglePeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	const int numPeaks = result.size();
	// There are always going to be two peaks (guaranteed by initialPeakGuess)
	
	//const bool nonZeroPedestal = (adcData[0] + adcData[1])*0.5 > 4.0 / TMath::Sqrt2() * 3.0;
	if (numPeaks == 1)
	{
		SinglePeakFit singlePeak(_initParams);
		singlePeak.process(adcData, initialGuess, result);
	}
	else // If there is a double peak
	{
		LXPeakFit doublePeak(_initParams);
		doublePeak.process(adcData, initialGuess, result);
	}
}
}
}
