// MAKE SURE THAT initParams is getting passed everywhere correctly 

#include "TrkChargeReco/inc/ComboPeakFit.hh"
#include "TrkChargeReco/inc/FitModelRoot.hh"


//SumADC::SumADC(const ConfigStruct &initParams) : PeakFitRootBase(initParams){}
namespace mu2e {

namespace TrkChargeReco {

void SumADC::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
 	int sum = 0;
	for (auto i = 0;  i < _initParams._numSamplesPerHit; ++i)
	{
		sum += adcData[i];
	}

	// Subtract the pedestal
	sum -= 2 * _initParams._defaultPedestal;
	
	// Unlike the other Find Peak Methods Sum ADC return the sum in result._peakHeight and does not alter the peak time 	
	
	resultantPeakData peakData;
	peakData._peakHeight = sum; 
	result[0] = peakData;
}

// SinglePeakFit normal constructor with ConfigStruct initilization parameters
SinglePeakFit::SinglePeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
  _fitModel = TF1("fitModel",FitModelRoot::singlePeakTrunc,0.0,_initParams._hitPeriod,3);
}

// Fills result using adc waveform data using by fitting with the early peak model
// NOTE : This function may begin with peak data provided in result which is replaced
void SinglePeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	// Set initial fit parameters
	const double timeshift = 30.0;
	const double scalingFactor = TMath::Max((initialGuess[0]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double sigma = 10.0;
	const Double_t initialParameters[3] = {timeshift, scalingFactor, sigma};

	Double_t finalParameters[3];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters); 
	fitParams2ResultantData(finalParameters, result);
}

void SinglePeakFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData peakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	peakData._peakTime = fitParameters[0];
	peakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	result[0] = peakData;
}	

// SinglePeakFloatingPedestalFit normal constructor with ConfigStruct initilization parameters
SinglePeakFloatingPedestalFit::SinglePeakFloatingPedestalFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
  _fitModel = TF1("fitModel",FitModelRoot::singlePeakFloatingPedestalTrunc,0.0,_initParams._hitPeriod,4);
}

// Fills result using adc waveform data using by fitting with the singlePeakFloatingPedestalTrunc model
void SinglePeakFloatingPedestalFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{	
	// Set initial fit parameters
	const double timeshift = 30.0;
	const double scalingFactor = TMath::Max((initialGuess[0]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double sigma = 10.0;
	const double verticalShift = (adcData[0] + adcData[1]) * 0.5;
	const Double_t initialParameters[4] = {timeshift, scalingFactor , verticalShift, sigma};

	Double_t finalParameters[4];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters); 
	fitParams2ResultantData(finalParameters, result);
}

void SinglePeakFloatingPedestalFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData peakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	peakData._peakTime = fitParameters[0];
	peakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	result[0] = peakData;
}	

// EXPeakFit normal constructor with ConfigStruct initilization parameters
EXPeakFit::EXPeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::EXPeakTrunc,0.0,_initParams._hitPeriod,4);
}

// Fills result using adc waveform data using by fitting with the EXPeakTrunc model
// NOTE : This function may begin with peak data provided in result which is replaced
void EXPeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	// Set initial fit parameters
	const double timeshift = 30.0;
	const double scalingFactor = TMath::Max((initialGuess[1]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double Q = initialGuess[0]._peakHeight;
	const double sigma = 10.0;
	Double_t initialParameters[4] = {timeshift, scalingFactor, Q, sigma};

	Double_t finalParameters[4];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void EXPeakFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{

	// Dynamic Pedestal is treated as the first peak with a "peak time" of 0.0
	resultantPeakData earlyPeakData;

	earlyPeakData._peakTime = 0.0;
	earlyPeakData._peakHeight = fitParameters[2];

	resultantPeakData peakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	peakData._peakTime = fitParameters[0];
	peakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	result[0] = earlyPeakData;
	result[1] = peakData;
}

LXPeakFit::LXPeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::LXPeakTrunc,0.0,_initParams._hitPeriod,4);
}

// Fills result using adc waveform data using by fitting with the LXPeakTrunc model
// NOTE : This function may begin with peak data provided in result which is replaced
void LXPeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	// Set initial fit parameters
	const double timeShift0 = 30.0;
	const double scalingFactor0 = TMath::Max((initialGuess[0]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double timeshift1 = initialGuess[1]._peakTime - initialGuess[0]._peakTime;
	const double scalingFactor1 = TMath::Max((initialGuess[1]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const Double_t initialParameters[5] = {timeShift0, scalingFactor0, timeshift1, scalingFactor1};

	Double_t finalParameters[5];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void LXPeakFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData firstPeakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	firstPeakData._peakTime = fitParameters[0];
	firstPeakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	resultantPeakData secondPeakData;
	secondPeakData._peakTime = fitParameters[2];
	secondPeakData._peakHeight = fitParameters[3] * _initParams._scalingFactor2bits;

	result[0] = firstPeakData;
	result[1] = secondPeakData;
}


LXPeakFloatingPedestalFit::LXPeakFloatingPedestalFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::LXPeakFloatingPedestalTrunc,0.0,_initParams._hitPeriod,5);
}

// Fills result using adc waveform data using by fitting with the LXPeakFloatingPedestalTrunc model
void LXPeakFloatingPedestalFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	// Set initial fit parameters
	const double timeShift0 = 30.0;
	const double scalingFactor0 = TMath::Max((initialGuess[0]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double verticalShift = (adcData[0] + adcData[1]) * 0.5; // This has implicit casting
	const double timeshift1 = initialGuess[1]._peakTime - initialGuess[0]._peakTime;
	const double scalingFactor1 = TMath::Max((initialGuess[1]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const Double_t initialParameters[5] = {timeShift0, scalingFactor0, verticalShift, timeshift1, scalingFactor1};

	Double_t finalParameters[5];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void LXPeakFloatingPedestalFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData firstPeakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	firstPeakData._peakTime = fitParameters[0];
	firstPeakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	resultantPeakData secondPeakData;
	secondPeakData._peakTime = fitParameters[3];
	secondPeakData._peakHeight = fitParameters[4] * _initParams._scalingFactor2bits;

	result[0] = firstPeakData;
	result[1] = secondPeakData;
}

// ELXPeakFit normal constructor with ConfigStruct initilization parameters
ELXPeakFit::ELXPeakFit(const ConfigStruct &initParams) : PeakFitRootBase(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::ELXPeakTrunc,0.0,_initParams._hitPeriod,5);
}

// Fills result using adc waveform data using by fitting with the ELXPeakTrunc model
// NOTE : This function may begin with peak data provided in result which is replaced
void ELXPeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	const double timeShift0 = 30.0;
	const double scalingFactor0 = TMath::Max((initialGuess[1]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);
	const double Q = initialGuess[0]._peakHeight;
	const double timeshift1 = initialGuess[2]._peakTime - initialGuess[1]._peakTime;
	const double scalingFactor1 = TMath::Max((initialGuess[2]._peakHeight - _initParams._defaultPedestal) * _initParams._bits2scalingFactor, 1000.0);

	Double_t initialParameters[5] = {timeShift0, scalingFactor0, Q, timeshift1, scalingFactor1};
	Double_t finalParameters[5];

	PeakFitRootBase::adcWaveform2TGraphErrors(adcData, _fitData);
	PeakFitRootBase::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void ELXPeakFit::fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result)
{
	// First peak is dynamic pedestal with a "peak time" of 0.0
	resultantPeakData earlyPeakData;

	earlyPeakData._peakTime = 0.0;
	earlyPeakData._peakHeight = fitParameters[2];


	resultantPeakData firstPeakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	firstPeakData._peakTime = fitParameters[0];
	firstPeakData._peakHeight = fitParameters[1] * _initParams._scalingFactor2bits;

	resultantPeakData secondPeakData;
	secondPeakData._peakTime = fitParameters[3];
	secondPeakData._peakHeight = fitParameters[4] * _initParams._scalingFactor2bits;

	result[0] = earlyPeakData;
	result[1] = firstPeakData;
	result[2] = secondPeakData;
}

// NOTE : This function may begin with peak data provided in result which is replaced
void ComboPeakFit::process(const adcWaveform adcData, const resultantHitData &initialGuess, resultantHitData &result)
{
	const int numPeaks = result.size();
	if (initialGuess[0]._peakTime == 0.0) // If there is a dynamic pedestal
	{	
		// There are always going to be two peaks (guaranteed by initialPeakGuess)
		if (numPeaks == 2)
		{
			EXPeakFit singlePeak(_initParams);
			singlePeak.process(adcData, initialGuess, result);
		}
		else if (numPeaks == 3)
		{
			ELXPeakFit LXPeak(_initParams);
			LXPeak.process(adcData, initialGuess, result);
		}
	}
	// If there is no dynamic pedestal
	// TODO : Finish implementing thie possibility of nonzero pedestal
	else
	{
		const bool nonZeroPedestal = (adcData[0] + adcData[1])*0.5 > 4.0 / TMath::Sqrt2() * 3.0;
		if (numPeaks == 1)
		{
			if (nonZeroPedestal)
			{
				SinglePeakFloatingPedestalFit singlePeak(_initParams);
				singlePeak.process(adcData, initialGuess, result);
			}
			else
			{
				SinglePeakFit singlePeak(_initParams);
				singlePeak.process(adcData, initialGuess, result);
			}
		}
		if (numPeaks == 2)
		{
			if (nonZeroPedestal)
			{
				LXPeakFloatingPedestalFit LXPeak(_initParams);
				LXPeak.process(adcData, initialGuess, result);
			}
			else		
			{
				LXPeakFit LXPeak(_initParams);
				LXPeak.process(adcData, initialGuess, result);
			}
		}
	}	
}
}
}


