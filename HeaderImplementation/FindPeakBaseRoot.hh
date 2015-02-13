#ifndef FindPeakBaseRoot_hh
#define FindPeakBaseRoot_hh

#include "FindPeakBase.hh"
#include "TF1.h"
#include "TGraphErrors.h"

class FindPeakBaseRoot : public FindPeakBase{
	public:

		// FindPeakBaseRoot normal constructor with ConfigStruct initilization parameters
		FindPeakBaseRoot(const ConfigStruct &initParams) : FindPeakBase(initParams){};

	protected:

		// Fits a model function to a waveform
		void fitModel2NormalizedWaveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters);

		// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
		void adcWaveform2TGraphErrors(const adcWaveform adcData, TGraphErrors &fitData);

		TGraphErrors _fitData;
		TF1 _fitModel;

};
#endif