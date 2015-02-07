// MAKE SURE THAT initParams is getting passed everywhere correctly 

#include "FindMultiplePeak.hh"
#include "FitModelRoot.hh"

// TODO : ADD CLASS FindSinglePeakWithConstantPedestal


// FindSinglePeak normal constructor with configStruct initilization parameters
FindSinglePeak::FindSinglePeak(const configStruct &initParams) : FindPeakBaseRoot(initParams)
{
  _fitModel = TF1("fitModel",FitModelRoot::convolutionSinglePeakWithConstantPedestal,0.0,_hitPeriod,4);
}


// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
// NOTE : This function may begin with peak data provided in result which is replaced
void FindSinglePeak::process(const adcWaveform adcData, resultantHitData &result)
{	
	// Set initial fit parameters
	const double timeshift = 30.0;
	const double scalingFactor = TMath::Max(result[0]._peakHeight * _bits2scalingFactor, 1000.0);
	std::cout << "peak height : " << result[0]._peakHeight << std::endl; 
	const double sigma = 10.0;
	double verticalShift = 0.0;
	const Double_t initialParameters[4] = {timeshift, scalingFactor , verticalShift, sigma};
	std::cout << "scalingFactor : " << scalingFactor << endl;

	// DEAL WITH CASE OF ZERO PEDESTAL
	//const bool nonZeroPedestal = (qAdc[0] + qAdc[1])*0.5 > 4.0 / sqrt(2.0) * 3.0;
	
	Double_t finalParameters[4];

	FindPeakBaseRoot::adcWaveform2TGraphErrors(adcData, _fitData);
	FindPeakBaseRoot::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	std::cout <<  " final param 0  : " << finalParameters[0] << std::endl; 
	std::cout <<  " final param 1  : " << finalParameters[1] << std::endl;
	std::cout <<  " final param 2  : " << finalParameters[2] << std::endl; 
	std::cout <<  " final param 3  : " << finalParameters[3] << std::endl;  
	fitParams2ResultantData(finalParameters, result);
}


void FindSinglePeak::fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData peakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	peakData._peakTime = fitParameters[0];
	peakData._peakHeight = fitParameters[1] * _scalingFactor2bits;

	result[0] = peakData;
}	


// FindSinglePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
FindSinglePeakWithDynamicPedestal::FindSinglePeakWithDynamicPedestal(const configStruct &initParams) : FindPeakBaseRoot(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::convolutionSinglePeakWithDynamicPedestal,0.0,_hitPeriod,4);
}

// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
// NOTE : This function may begin with peak data provided in result which is replaced
void FindSinglePeakWithDynamicPedestal::process(const adcWaveform adcData, resultantHitData &result)
{

// Set initial fit parameters
const double timeshift = 30.0;
const double scalingFactor = result[1]._peakHeight * _bits2scalingFactor;
const double Q = result[0]._peakHeight;
const double sigma = 10.0;
Double_t initialParameters[4] = {timeshift, scalingFactor, Q, sigma};

Double_t finalParameters[4];

FindPeakBaseRoot::adcWaveform2TGraphErrors(adcData, _fitData);
FindPeakBaseRoot::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
fitParams2ResultantData(finalParameters, result);
}

void FindSinglePeakWithDynamicPedestal::fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
{

	// Dynamic Pedestal is treated as the first peak with a "peak time" of 0.0
	resultantPeakData dynamicPedestalData;

	dynamicPedestalData._peakTime = 0.0;
	dynamicPedestalData._peakHeight = fitParameters[2];

	resultantPeakData peakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	peakData._peakTime = fitParameters[0];
	peakData._peakHeight = fitParameters[1] *_scalingFactor2bits;

	result[0] = dynamicPedestalData;
	result[1] = peakData;

}


FindDoublePeak::FindDoublePeak(const configStruct &initParams) : FindPeakBaseRoot(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::doublePeak,0.0,_hitPeriod,5);
}

// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
// NOTE : This function may begin with peak data provided in result which is replaced
void FindDoublePeak::process(const adcWaveform adcData, resultantHitData &result)
{
	// Set initial fit parameters
	const double timeShift0 = 30.0;
	const double scalingFactor0 = TMath::Max(result[0]._peakHeight / 0.015, 1000.0);
	double verticalShift = (adcData[0] + adcData[1]) * 0.5; // This has implicit casting
	const double timeshift1 = result[1]._peakTime - result[0]._peakTime;
	const double scalingFactor1 = TMath::Max(result[1]._peakHeight / 0.015, 1000.0);
	const Double_t initialParameters[5] = {timeShift0, scalingFactor0, verticalShift, timeshift1, scalingFactor1};

	Double_t finalParameters[5];

	FindPeakBaseRoot::adcWaveform2TGraphErrors(adcData, _fitData);
	FindPeakBaseRoot::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void FindDoublePeak::fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
{
	resultantPeakData firstPeakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	firstPeakData._peakTime = fitParameters[0];
	firstPeakData._peakHeight = fitParameters[1] * _scalingFactor2bits;

	resultantPeakData secondPeakData;
	secondPeakData._peakTime = fitParameters[3];
	secondPeakData._peakHeight = fitParameters[4] * _scalingFactor2bits;

	result[0] = firstPeakData;
	result[1] = secondPeakData;
}



// FindDoublePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
FindDoublePeakWithDynamicPedestal::FindDoublePeakWithDynamicPedestal(const configStruct &initParams) : FindPeakBaseRoot(initParams)
{
	_fitModel = TF1("fitModel",FitModelRoot::doublePeakWithDynamicPedestal,0.0,_hitPeriod,5);
}

// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
// NOTE : This function may begin with peak data provided in result which is replaced
void FindDoublePeakWithDynamicPedestal::process(const adcWaveform adcData, resultantHitData &result)
{
	const double timeShift0 = 30.0;
	const double scalingFactor0 = TMath::Max(result[1]._peakHeight * _bits2scalingFactor, 1000.0);
	const double Q = result[0]._peakHeight;
	const double timeshift1 = result[2]._peakTime - result[1]._peakTime;
	const double scalingFactor1 = TMath::Max(result[2]._peakHeight * _bits2scalingFactor, 1000.0);

	Double_t initialParameters[5] = {timeShift0, scalingFactor0, Q, timeshift1, scalingFactor1};
	Double_t finalParameters[5];

	FindPeakBaseRoot::adcWaveform2TGraphErrors(adcData, _fitData);
	FindPeakBaseRoot::fitModel2NormalizedWaveform(_fitModel, _fitData, initialParameters, finalParameters);
	fitParams2ResultantData(finalParameters, result);
}

void FindDoublePeakWithDynamicPedestal::fitParams2ResultantData(Double_t *fitParameters, resultantHitData &result)
{
	// First peak is dynamic pedestal with a "peak time" of 0.0
	resultantPeakData dynamicPedestalData;

	dynamicPedestalData._peakTime = 0.0;
	dynamicPedestalData._peakHeight = fitParameters[2];


	resultantPeakData firstPeakData;

	// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
	firstPeakData._peakTime = fitParameters[0];
	firstPeakData._peakHeight = fitParameters[1] * _scalingFactor2bits;

	resultantPeakData secondPeakData;
	secondPeakData._peakTime = fitParameters[3];
	secondPeakData._peakHeight = fitParameters[4] * _scalingFactor2bits;

	result[0] = dynamicPedestalData;
	result[1] = firstPeakData;
	result[2] = secondPeakData;
}

// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
// NOTE : This function may begin with peak data provided in result which is replaced
void FindMultiplePeaks::process(const adcWaveform adcData, resultantHitData &result)
{
	FindPeakBaseRoot::adcWaveform2TGraphErrors(adcData, _fitData);
	findPeaks(_fitData, initParams, result);
	std::cout << "Peak Time : " << result[0]._peakTime << std::endl;
	std::cout << "Peak Height : " << result[0]._peakHeight << std::endl;


	const int numPeaks = result.size();

	if (result[0]._peakTime == 0.0) // If there is a dynamic pedestal
	{
		// If the only peak is a dynamic pedestal 
		// search for another peak
		if (numPeaks == 1) 
		{
			dynamicPedestalAddPeak(_fitData, result);
			FindSinglePeakWithDynamicPedestal singlePeak(initParams);
			singlePeak.process(adcData, result);
		}
		else if (numPeaks == 2)
		{
			FindSinglePeakWithDynamicPedestal singlePeak(initParams);
			singlePeak.process(adcData, result);
		}
		else if (numPeaks == 3)
		{
			FindDoublePeakWithDynamicPedestal doublePeak(initParams);
			doublePeak.process(adcData, result);
		}
	}
	// If there is no dynamic pedestal
	else
	{
		if (numPeaks == 1)
		{
			FindSinglePeak singlePeak(initParams);
			singlePeak.process(adcData, result);
		}
		if (numPeaks == 2)
		{
			FindDoublePeak doublePeak(initParams);
			doublePeak.process(adcData, result);
		}

	}	
}

// TODO : FIGURE OUT WHY THIS FUNCTION IS NO LONGER RETURNING THE CORRECT VALUE
// Performs explicit peak search on adc waveform data
void FindMultiplePeaks::findPeaks(TGraphErrors &gr, const configStruct &initParams, resultantHitData &result, double sigma)
{
		int ientry = 0; // Start time at 0
		const double *measurementTimes = gr.GetX();
		const double *adcValues = gr.GetY();

		while(ientry < initParams._numSamplesPerHit)
		{
		double adcValue = adcValues[ientry];
		double tMax = measurementTimes[ientry];
		double adcMax = adcValue;
		double adcPrev = adcValue;

		int jentry = ientry + 1;
		bool descending = false;
		while (jentry < initParams._numSamplesPerHit)
		{
			adcValue = adcValues[jentry];
			descending |= ((adcPrev-adcValue) > (TMath::Sqrt2()*initParams._adcError*sigma));

			if (descending && (adcValue-adcPrev > (TMath::Sqrt2()*initParams._adcError*sigma)))
			{
			break;
			}
			else
			{
			if (adcValue > adcMax)
			{
  			adcMax  = adcValue;
  			tMax = measurementTimes[jentry];
				}
			adcPrev = adcValue;
			ientry = jentry;
			++jentry;
			}
			}
		resultantPeakData peakData(tMax, adcMax);
		result.push_back(peakData);
		++ientry;
		}
}

// This function searches for another peak in the waveform data by subtracting out a dynamic pedestal 
// from the adc waveform and finding the maximum adc value in the "subtracted data".
// This function is applied when no peak is found in the explicit peak search (findPeaks).
void FindMultiplePeaks::dynamicPedestalAddPeak(TGraphErrors &gr, resultantHitData &result)
{	
	// This maybe could be done using linear algebra vectors
	// instead of arrays
	const Double_t *adcValues = gr.GetY();
	const Double_t *measurementTimes = gr.GetX();
	Double_t subtractedValues[FindPeakBase::_initParams._numSamplesPerHit];

	Double_t dynamicPedstalParam[1] = {adcValues[0]};
	Double_t dynamicPedestalX[1];

	for (int i = 0; i < FindPeakBase::_initParams._numSamplesPerHit; ++i)
	{
		dynamicPedestalX[0] = measurementTimes[i];
		subtractedValues[i] = adcValues[i] - dynamicPedestal(dynamicPedestalX, dynamicPedstalParam);
	}

	// New peak is max value of difference between of adc values and dynamic pedestal
	const Float_t newAdcPeak = TMath::MaxElement(FindPeakBase::_initParams._numSamplesPerHit, subtractedValues);
	const Float_t newTPeak = TMath::LocMax(FindPeakBase::_initParams._numSamplesPerHit, subtractedValues);

	resultantPeakData newPeakData(newTPeak, newAdcPeak);
	result.push_back(newPeakData);
}

