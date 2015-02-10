#ifndef FindPeakBaseRoot_hh
#define FindPeakBaseRoot_hh

#include "FindPeakBase.hh"

class FindPeakBaseRoot : public FindPeakBase{
	public:
		
		// Fills result using adc waveform data
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result) = 0;

		// FindPeakBaseRoot normal constructor with configStruct initilization parameters
		FindPeakBaseRoot(const configStruct &initParams) : FindPeakBase(initParams){};

	protected:

		// Fits a model function to a waveform
		void fitModel2NormalizedWaveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters);

		// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
		void adcWaveform2TGraphErrors(adcWaveform adcData, TGraphErrors &fitData);

};
#endif