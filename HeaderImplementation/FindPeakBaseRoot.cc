#include "FindPeakBaseRoot.hh"

// Fits a model function to a waveform
// Note: It is assumed that the waveform being fitted is normalized. 
// Thus, the scaling factor parameter which is passed in must be not in units of bits but scaled to this normaliztion
// This can be done by multiplying by the _bits2scalingfactor conversion factor
void FindPeakBaseRoot::fitModel2NormalizedWaveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters)
{
	// These lines will be replaced with the chi-square minimization
	TF1 *fitModelPtr = &fitModel; 
	TGraphErrors *fitDataPtr = &fitData;
	fitModel.SetParameters(initialParameters);
	fitDataPtr->Fit(fitModelPtr,"QN");

	const Int_t numParameters = fitModel.GetNumberFreeParameters();

	for (int i = 0; i < numParameters; ++i)
	{
		fitParameters[i] = fitModel.GetParameter(i);
		std::cout << "init param [" << i << "]" << initialParameters[i] << std::endl;
	}
}


// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
void FindPeakBaseRoot::adcWaveform2TGraphErrors(const adcWaveform adcData, TGraphErrors &fitData)
{
	Double_t adcDataTemp[_initParams._numSamplesPerHit];
	Double_t measurementTimes[_initParams._numSamplesPerHit];
	Double_t measurementTimesErrors[_initParams._numSamplesPerHit];
	Double_t adcDataErrors[_initParams._numSamplesPerHit];

	for (int i = 0; i < _initParams._numSamplesPerHit; ++i)
	{
		adcDataTemp[i] = (Double_t) adcData[i];
		measurementTimes[i] = (Double_t) i * _initParams._measurementFrequency; 
		measurementTimesErrors[i] = 0.0;
		adcDataErrors[i] = _initParams._adcError;
	}

	fitData = TGraphErrors(_initParams._numSamplesPerHit,measurementTimes,adcDataTemp,measurementTimesErrors,adcDataErrors);
}
