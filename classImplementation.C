#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>
#include "FitModel.C"

typedef unsigned int * adcWaveform;
typedef std::vector<Double_t*> resultantData;

struct paramStruct{
		const Double_t shapingTime;
		std::vector<unsigned int> *adcPeaks;
		std::vector<unsigned int> *adcPeakTimes;
		const Int_t numSamplesPerHit;
		const Double_t *adcErrors;
		const Double_t *timeMeasurementErrors;
		const Double_t *measurementTimes;
		const Double_t measurementFrequency; // ns 

		paramStruct() : shapingTime(25.0), adcPeaks(0),adcPeakTimes(0),numSamplesPerHit(8), measurementFrequency(20.0)
		{

			// Initialize ADC Errors and time measurement errors
			Double_t adcErrorsTemp[numSamplesPerHit];
			Double_t timeMeasurementErrorsTemp[numSamplesPerHit];
			Double_t measurementTimesTemp[numSamplesPerHit];

			unsigned int adcDataTemp[numSamplesPerHit];
			for (int i = 0; i < numSamplesPerHit; ++i)
			{
				adcErrorsTemp[i] = 3.0;
				timeMeasurementErrorsTemp[i] = 0.0;
				measurementTimesTemp[i] = i * measurementFrequency;
			}	
			adcErrors = adcErrorsTemp;
			timeMeasurementErrors = timeMeasurementErrorsTemp;
			measurementTimes = measurementTimesTemp;
		}

	};


class FindPeakBase{
	public:
		//virtual resultantData process(adcWaveform waveformData, std::vector<unsigned int> *adcPeaks, std::vector<unsigned int> *adcPeakTimes) = 0;
		virtual resultantData process() = 0;

		virtual ~FindPeakBase(){}

		FindPeakBase(){}

		FindPeakBase(paramStruct initParams, adcWaveform adcData)
		{
		// These 3 lines NEED to be deleted. 
		// They are a result of taking adc data from initParams instead of adcData
		Double_t adcDataTemp[8];

		for (int i = 0; i < 8; ++i)
		{
			adcDataTemp[i] = (Double_t) adcData[i];
		}


		fitData = new TGraphErrors(initParams.numSamplesPerHit,
									initParams.measurementTimes,
									adcDataTemp,
									initParams.timeMeasurementErrors,
								   	initParams.adcErrors);

		//this->initParams = initParams;
		this->adcData = adcData;
		}

	protected:
		TGraphErrors *fitData;
		TF1 *fitModel;
		paramStruct initParams; 
		adcWaveform adcData;

		// These should probably change from Float_t to Int_t (or unsigned int)
		std::vector<Float_t> adcPeaks;
		std::vector<Float_t> timePeaks;

		void fitModel2Data(TF1 *fitModel, TGraph *fitData, Double_t *initialParameters, Double_t *fitParameters)
		{
			// These lines will be replaced with the chi-square minimization
			fitModel->SetParameters(fitParameters);
			fitData->Fit(fitModel,"QN");
			fitParameters = fitModel->GetParameters();
		}
};

class FindSinglePeakWithDynamicPedestal : public FindPeakBase{
	public:
		FindSinglePeakWithDynamicPedestal(paramStruct initParams, adcWaveform adcData) : FindPeakBase(initParams, adcData){}

		virtual resultantData process()
		{
			convolutionSinglePeakWithDynamicPedestal *c = new convolutionSinglePeakWithDynamicPedestal();
			const double timeshift = 30.0;
			const double scalingFactor = (double) adcPeaks[1] / 0.015;
			const double Q = (double) adcPeaks[0];
			const double sigma = 10.0;

			fitModel = new TF1("fitModel",c,&convolutionSinglePeakWithDynamicPedestal::fitModel
								,0.0,140.0,4,"convolutionSinglePeakWithDynamicPedestal","fitModel");
			Double_t initialParameters[4] = {timeshift, scalingFactor, Q, sigma};
			Double_t finalParameters[4];

			fitModel2Data(fitModel, fitData, initialParameters, finalParameters);
			resultantData result;

			Double_t peakData[2] = {finalParameters[0], finalParameters[1]};
			result.push_back(peakData);
			return result;
		}
};

class FindSinglePeak : public FindPeakBase{
	public:
		FindSinglePeak(paramStruct initParams, adcWaveform adcData) : FindPeakBase(initParams, adcData){}

		virtual resultantData process()
		{
			convolutionSinglePeak *c = new convolutionSinglePeak();
			const double timeshift = 30.0;
			const double scalingFactor = TMath::Max(adcPeaks[0] / 0.015, 1000.0);
			const double sigma = 10.0;

			double verticalShift = 0.0;

			// DEAL WITH CASE OF ZERO PEDESTAL
			// const bool nonZeroPedestal = (qAdc[0] + qAdc[1])*0.5 > 4.0 / sqrt(2.0) * 3.0;

			fitModel = new TF1("fitModel",c,&convolutionSinglePeak::fitModel,0.0,140.0,4,"convolutionSinglePeak","fitModel");
			Double_t initialParameters[4] = {timeshift, scalingFactor, verticalShift, sigma};
			Double_t finalParameters[4];

			fitModel2Data(fitModel, fitData, initialParameters, finalParameters);
			resultantData result;

			Double_t peakData[2] = {finalParameters[0], finalParameters[1]};
			result.push_back(peakData);
			return result;
		}
};

class FindDoublePeak : public FindPeakBase{
	public:
		FindDoublePeak(paramStruct initParams, adcWaveform adcData) : FindPeakBase(initParams, adcData){}

		virtual resultantData process()
		{
			doublePeak *c = new doublePeak();
			fitModel = new TF1("fitModel",c,&doublePeak::fitModel,0.0,140.0,5,"doublePeak","fitModel");

			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(adcPeaks[0] / 0.015, 1000.0);
			double verticalShift = (adcData[0] + adcData[1]) * 0.5; // This has implicit casting
			const double timeshift1 = timePeaks[1] - timePeaks[0];
			const double scalingFactor1 = TMath::Max(adcPeaks[1] / 0.015, 1000.0);

			Double_t initialParameters[5] = {timeShift0, scalingFactor0, verticalShift, timeshift1, scalingFactor1};
			Double_t finalParameters[5];

			fitModel2Data(fitModel, fitData, initialParameters, finalParameters);
			resultantData result; 
			fitParams2ResultantData(result, finalParameters);
			return result;

		}

	protected:

		void fitParams2ResultantData(resultantData &result, Double_t *fitParameters)
		{
			Double_t firstPeakData[2];

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			firstPeakData[0] = fitParameters[0];
			firstPeakData[1] = fitParameters[1];

			Double_t secondPeakData[2];
			secondPeakData[0] = fitParameters[3];
			secondPeakData[1] = fitParameters[4];

			result.push_back(firstPeakData);
			result.push_back(secondPeakData);
		}
};


class FindDoublePeakWithDynamicPedestal : public FindDoublePeak{
	public:
		FindDoublePeakWithDynamicPedestal(paramStruct initParams, adcWaveform adcData) : FindDoublePeak(initParams, adcData){}

		virtual resultantData process()
		{
			doublePeakWithDynamicPedestal *c = new doublePeakWithDynamicPedestal();
			fitModel = new TF1("fitModel",c,&doublePeak::fitModel,0.0,140.0,5,"doublePeak","fitModel");
			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(adcPeaks[1] / 0.015, 1000.0);
			const double Q = adcPeaks[0];
			const double timeshift1 = timePeaks[2] - timePeaks[1];
			const double scalingFactor1 = TMath::Max(adcPeaks[2] / 0.015, 1000.0);

			Double_t initialParameters[5] = {timeShift0, scalingFactor0, Q, timeshift1, scalingFactor1};
			Double_t finalParameters[5];

			fitModel2Data(fitModel, fitData, initialParameters, finalParameters);
			resultantData result; 
			FindDoublePeak::fitParams2ResultantData(result, finalParameters);
			return result;
		}

};




class FindMultiplePeaks : public FindPeakBase{
	public:

		FindMultiplePeaks(paramStruct initParams, adcWaveform adcData)
		: FindPeakBase(initParams, adcData)
		{
			findPeaks(fitData,adcPeaks,timePeaks);
		}

		virtual resultantData process()
		{

			resultantData result;

			if (timePeaks[0] == 0.0) // If there is a dynamic pedestal
			{

				// If the only peak is a dynamic pedestal 
				// search for another peak
				if (timePeaks.size() == 1) 
				{
					dynamicPedestalAddPeak(fitData, timePeaks, adcPeaks);
					FindSinglePeakWithDynamicPedestal *singlePeak = new FindSinglePeakWithDynamicPedestal(initParams, adcData);
					result = singlePeak->process();
				}
				else if (timePeaks.size() == 2)
				{
					FindSinglePeakWithDynamicPedestal *singlePeak = new FindSinglePeakWithDynamicPedestal(initParams, adcData);
					result = singlePeak->process();
				}
				else if (timePeaks.size() == 3)
				{
					FindDoublePeakWithDynamicPedestal *doublePeak = new FindDoublePeakWithDynamicPedestal(initParams, adcData);
					result = doublePeak->process();
				}
			}
			// If there is no dynamic pedestal
			else
			{
				if (timePeaks.size() == 1)
				{
					FindSinglePeak *singlePeak = new FindSinglePeak(initParams, adcData);
					result = singlePeak->process();
				}
				if (timePeaks.size() == 2)
				{
					FindDoublePeak *doublePeak = new FindDoublePeak(initParams, adcData);
					result = doublePeak->process();
				}

			}	
			return result;
		}


		//resultantData process(adcWaveform waveformData, std::vector<unsigned int> *adcPeaks, std::vector<unsigned int> *adcPeakTimes)

		void findPeaks(TGraphErrors *gr,std::vector<Float_t>& tPeak, std::vector<Float_t>& adcPeak, double sigma = 3.0)
		{
  			int ientry = 0; // Start time at 0
  			const int nEntries = gr->GetN();
  			const double *measurementTimes = gr->GetX();
  			const double *adcValues = gr->GetY();
  			const double measurementError = gr->GetErrorY(0);

  			while(ientry < nEntries)
  			{
    			double adcValue = adcValues[ientry];
				double tMax = measurementTimes[ientry];
    			double adcMax = adcValue;
    			double adcPrev = adcValue;

    			int jentry = ientry + 1;
    			bool descending = false;
    			while (jentry < nEntries)
				{
					adcValue = adcValues[jentry];
      				descending |= ((adcPrev-adcValue) > (TMath::Sqrt(2.0)*measurementError*sigma));

      				if (descending && (adcValue-adcPrev > (TMath::Sqrt(2.0)*measurementError*sigma)))
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
      			tPeak.push_back(tMax);
      			adcPeak.push_back(adcMax);
      			++ientry;
  			}
		}

	private:
		void dynamicPedestalAddPeak(TGraphErrors *gr, std::vector<Float_t>& tPeak, std::vector<Float_t>& adcPeak)
		{	
			const int numSamples = gr->GetN();

			// This maybe could be done using linear algebra vectors
			// instead of arrays
			const Double_t *adcValues = gr->GetY();
			const Double_t *measurementTimes = gr->GetX();
			Double_t subtractedValues[numSamples];

			const Double_t dynamicPedstalValue = adcValues[0];

			// GET RID OF THIS
			const double shapingTime = 25.0; 

			for (int i = 0; i < gr->GetN(); ++i)
			{
				subtractedValues[i] = adcValues[i] - dynamicPedstalValue * exp(-measurementTimes[i] / shapingTime);
			}
			const Float_t newAdcPeak = TMath::MaxElement(numSamples,subtractedValues);
			const Float_t newTPeak = TMath::LocMax(numSamples, subtractedValues);

			adcPeak.push_back(newAdcPeak);
			tPeak.push_back(newTPeak);
		}
};

