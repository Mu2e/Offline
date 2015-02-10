#ifndef FindMultiplePeak_hh
#define FindMultiplePeak_hh

#include "FindPeakBaseRoot.hh"

class FindSinglePeak : public FindPeakBaseRoot{
	public:
		// FindSinglePeak normal constructor with configStruct initilization parameters
		FindSinglePeak(const configStruct &initParams);

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result);

	private:
		void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);

		TGraphErrors _fitData;
		TF1 _fitModel;
};

class FindSinglePeakWithDynamicPedestal : public FindPeakBaseRoot{
	public:

		// FindSinglePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindSinglePeakWithDynamicPedestal(const configStruct &initParams);

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result);

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;

		void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);

};

class FindDoublePeak : public FindPeakBaseRoot{
	public:
		FindDoublePeak(const configStruct &initParams);

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result);

	protected:

		void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;
};

class FindDoublePeakWithDynamicPedestal : public FindPeakBaseRoot{
	public:
		// FindDoublePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindDoublePeakWithDynamicPedestal(const configStruct &initParams);

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result);

	private:
		TGraphErrors _fitData;
		TF1 _fitModel;


		void fitParams2ResultantData(const Double_t *fitParameters, resultantHitData &result);
};

class FindMultiplePeaks : public FindPeakBaseRoot{
	public:

		// FindMultiplePeaks normal constructor with configStruct initilization parameters
		FindMultiplePeaks(const configStruct &initParams) : FindPeakBaseRoot(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(const adcWaveform adcData, resultantHitData &result);

	private:
		// Performs explicit peak search on adc waveform data
		void findPeaks(const TGraphErrors &gr, const configStruct &initParams, resultantHitData &result, const double sigma = 3.0);

		// This function searches for another peak in the waveform data by subtracting out a dynamic pedestal 
		// from the adc waveform and finding the maximum adc value in the "subtracted data".
		// This function is applied when no peak is found in the explicit peak search (findPeaks).
		void dynamicPedestalAddPeak(const TGraphErrors &gr, resultantHitData &result);

		TGraphErrors _fitData;
		TF1 _fitModel;


};
#endif
