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
	 	adcWaveform adcData;

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

				// Dummy value	
				adcDataTemp[i] = 2;
			}	
			adcErrors = adcErrorsTemp;
			timeMeasurementErrors = timeMeasurementErrorsTemp;
			measurementTimes = measurementTimesTemp;
			adcData = adcDataTemp;
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
		Double_t adcTemp[8];
		for (int i = 0; i < 8; ++i)
		{adcTemp[i] = initParams.adcData[i];}

		fitData = new TGraphErrors(initParams.numSamplesPerHit,
									initParams.measurementTimes,
									adcTemp,
									initParams.timeMeasurementErrors,
								   	initParams.adcErrors);
		}

	protected:
		TGraphErrors *fitData;
		TF1 *fitModel;

		// These should probably change from Float_t to Int_t (or unsigned int)
		std::vector<Float_t> adcPeaks;
		std::vector<Float_t> timePeaks;

		void fitModel2Data(TF1 *fitModel, TGraph *fitData, Double_t *fitParameters)
		{
			fitData->Fit(fitModel,"QN");
			fitParameters = fitModel->GetParameters();
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
			if (timePeaks[0] == 0.0) // If there is a dynamic pedestal
			{

				// If the only peak is a dynamic pedestal 
				// search for another peak
				if (timePeaks.size() == 1) 
				{
					dynamicPedestalAddPeak(fitData, timePeaks, adcPeaks);
					FindSinglePeakWithDynamicPedestal *singlePeak = new FindSinglePeakWithDynamicPedestal(initParams, adcData);
					singlePeak->process();
				}
				else if (timePeaks.size() == 2)
				{
					FindSinglePeakWithDynamicPedestal *singlePeak = new FindSinglePeakWithDynamicPedestal(initParams, adcData);
					singlePeak->process();
				}
				else if (timePeaks.size() == 3)
				{
					FindDoublePeakWithDynamicPedestal *doublePeak = new FindDoublePeakWithDynamicPedestal(initParams, adcData);
					doublePeak->process();
				}
			}
			// If there is no dynamic pedestal
			else
			{
				if (timePeaks.size() == 1)
				{
					FindSinglePeak *singlePeak = new FindSinglePeak(initParams, adcData);
				}
				if (timePeaks.size() == 2)
				{
					FindDoublePeak *doublePeak = new FindDoublePeak(initParams, adcData);
				}

			}

			
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

class FindSinglePeak : public FindPeakBase{
	public:
		/**virtual resultantData process(adcWaveform waveformData,std::vector<unsigned int> *adcPeaks, std::vector<unsigned int> *adcPeakTimes)
		{
			dummy[0] = 1.0;
			dummy[1] = 2.0;
			dummy[2] = 3.0;
			resultantData returningVector;
			returningVector.push_back(dummy);
			return returningVector;
		}**/
};

class FindSinglePeakWithDynamicPedestal : public FindPeakBase{
	public:
};

class FindDoublePeakWithDynamicPedestal : public FindPeakBase{
	public:
};

class FindDoublePeak : public FindPeakBase{
};

