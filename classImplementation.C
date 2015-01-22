#include <vector>
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>
#include "FitModel.C"

typedef unsigned int * adcWaveform;


struct resultantPeakData
{
	Float_t peakHeight; // in units of bits
	Float_t peakTime; // time of peak relative to 140.0 ns interval

	// Default constructor - should probably be deleted
	resultantPeakData() : peakTime(0.0), peakHeight(0.0){};

	// Check the order of the parameters where this is called
	resultantPeakData(Float_t peakTime, Float_t peakHeight) : peakTime(peakTime), peakHeight(peakHeight){};
};

// This struct contains all parameters which remain constant throughout the simulation
struct configStruct{
    const Double_t shapingTime; // Shaping time (in units of ns)
    const Int_t numSamplesPerHit; // Number of samples measured per hit
    const Double_t adcError; // Assumes constant error for all adc measurements (in units of bits)
    const Double_t measurementFrequency; // Sample frequency of adc values (in units of ns)
    const Double_t truncationLevel; // Level of truncation of waveform (in units of bits)
    const Double_t defaultPedestal; // Count value corresponding to the default pedestal (in units of bits)

    configStruct() : shapingTime(25.0), 
             numSamplesPerHit(8), 
             measurementFrequency(20.0), 
             adcError(3.0), 
             truncationLevel(1023.0),
             defaultPedestal(64.0){}
}; 

// This is object top which will be filled by the process method 
typedef std::vector<resultantPeakData> resultantHitData;


// Virtual class providing structure for FindSinglePeak, FindDoublePeak, FindMutiplePeaks, etc. 
class FindPeakBase{
	public:
		
		// Fills result using adc waveform data
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData) = 0;

		// Destructor
		virtual ~FindPeakBase(){}

		// Default Constructor
		FindPeakBase(){}

		// FindPeakBase normal constructor with configStruct initilization parameters
		FindPeakBase(const configStruct &initParams) : initParams(initParams){}

	protected:

		const configStruct initParams; 

		// These should probably change from Float_t to Int_t (or unsigned int)


		// Fits a model function to a waveform
		void fitModel2Waveform(TF1 &fitModel, TGraphErrors &fitData, const Double_t *initialParameters, Double_t *fitParameters)
		{
 			// These lines will be replaced with the chi-square minimization
			TF1 * fitModelPtr = &fitModel; 
			TGraphErrors *fitDataPtr = &fitData;
			fitModel.SetParameters(initialParameters);
			fitDataPtr->Fit(fitModelPtr,"QN");

			const Int_t numParameters = fitModel.GetNumberFreeParameters();
			std::cout << "numParam : " << numParameters << std::endl;
			for (int i = 0; i < numParameters; ++i)
			{
				fitParameters[i] = fitModel.GetParameter(i);
			}
		}

		// Converts adcWaveform object to TGraphErrors object for easier manipulation in ROOT
		void adcWaveform2TGraphErrors(adcWaveform adcData, TGraphErrors &fitData)
		{
			Double_t adcDataTemp[initParams.numSamplesPerHit];
			Double_t measurementTimes[initParams.numSamplesPerHit];
			Double_t measurementTimesErrors[initParams.numSamplesPerHit];
			Double_t adcDataErrors[initParams.numSamplesPerHit];

			for (int i = 0; i < initParams.numSamplesPerHit; ++i)
			{
				adcDataTemp[i] = (Double_t) adcData[i];
				measurementTimes[i] = (Double_t) i * initParams.measurementFrequency; 
				measurementTimesErrors[i] = 0.0;
				adcDataErrors[i] = initParams.adcError;
			}

			fitData = TGraphErrors(initParams.numSamplesPerHit,adcDataTemp,measurementTimes,measurementTimesErrors,adcDataErrors);
		}

		// Precomputed constants
		const Double_t bits2scalingFactor = initParams.shapingTime * TMath::E(); // approximately 67.96
		const Double_t scalingFactor2bits  = 1.0 / bits2scalingFactor; // approximately 0.0147
		const Double_t hitPeriod = initParams.numSamplesPerHit * (initParams.measurementFrequency - 1.0);

};

class FindSinglePeakWithDynamicPedestal : public FindPeakBase{
	public:

		// FindSinglePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindSinglePeakWithDynamicPedestal(const configStruct &initParams) : FindPeakBase(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData)
		{

			// Set initial fit parameters
			const double timeshift = 30.0;
			const double scalingFactor = result[1].peakHeight * bits2scalingFactor;
			const double Q = result[0].peakHeight;
			const double sigma = 10.0;
			Double_t initialParameters[4] = {timeshift, scalingFactor, Q, sigma};

			// Define fit function
			fitModel = TF1("fitModel",convolutionSinglePeakWithDynamicPedestal,0.0,hitPeriod,4);

			Double_t finalParameters[4];

			adcWaveform2TGraphErrors(adcData, fitData);
			fitModel2Waveform(fitModel, fitData, initialParameters, finalParameters);
			fitParams2ResultantData(result,finalParameters);
		}

	private:
		TGraphErrors fitData;
		TF1 fitModel;

		void fitParams2ResultantData(resultantHitData &result, Double_t *fitParameters)
		{
			resultantPeakData peakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			peakData.peakTime = fitParameters[0];
			peakData.peakHeight = fitParameters[1];

			result[0] = peakData;

			// Since dynamic pedestal is not counted as peak
			result.pop_back();
		}




};

class FindSinglePeak : public FindPeakBase{
	public:
		// FindSinglePeak normal constructor with configStruct initilization parameters
		FindSinglePeak(const configStruct &initParams) : FindPeakBase(initParams){}


		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData)
		{	
			// Set initial fit parameters
			const double timeshift = 30.0;
			const double scalingFactor = TMath::Max(result[0].peakHeight * bits2scalingFactor, 1000.0);
			const double sigma = 10.0;
			double verticalShift = 0.0;
			const Double_t initialParameters[4] = {timeshift, scalingFactor, verticalShift, sigma};

			
			// DEAL WITH CASE OF ZERO PEDESTAL
			//const bool nonZeroPedestal = (qAdc[0] + qAdc[1])*0.5 > 4.0 / sqrt(2.0) * 3.0;

			// Define fit function
			fitModel = TF1("fitModel",convolutionSinglePeakWithConstantPedestal,0.0,hitPeriod,4);
	
			Double_t finalParameters[4];

			adcWaveform2TGraphErrors(adcData, fitData);
			fitModel2Waveform(fitModel, fitData, initialParameters, finalParameters);
			fitParams2ResultantData(result,finalParameters);
		}
	private:
		void fitParams2ResultantData(resultantHitData &result, Double_t *fitParameters)
		{
			resultantPeakData peakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			peakData.peakTime = fitParameters[0];
			peakData.peakHeight = fitParameters[1];

			result[0] = peakData;
		}



		TGraphErrors fitData;
		TF1 fitModel;
};

class FindDoublePeak : public FindPeakBase{
	public:
		FindDoublePeak(const configStruct &initParams) : FindPeakBase(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData)
		{
			// Set initial fit parameters
			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(result[0].peakHeight / 0.015, 1000.0);
			double verticalShift = (adcData[0] + adcData[1]) * 0.5; // This has implicit casting
			const double timeshift1 = result[1].peakTime - result[0].peakTime;
			const double scalingFactor1 = TMath::Max(result[1].peakHeight / 0.015, 1000.0);
			const Double_t initialParameters[5] = {timeShift0, scalingFactor0, verticalShift, timeshift1, scalingFactor1};

			// Define fit function
			fitModel = TF1("fitModel",doublePeak,0.0,hitPeriod,5);

			Double_t finalParameters[5];

			adcWaveform2TGraphErrors(adcData, fitData);
			fitModel2Waveform(fitModel, fitData, initialParameters, finalParameters);
			fitParams2ResultantData(result, finalParameters);
		}

	protected:

		void fitParams2ResultantData(resultantHitData &result, Double_t *fitParameters)
		{
			resultantPeakData firstPeakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			firstPeakData.peakTime = fitParameters[0];
			firstPeakData.peakHeight = fitParameters[1];

			resultantPeakData secondPeakData;
			secondPeakData.peakTime = fitParameters[3];
			secondPeakData.peakHeight = fitParameters[4];

			result[0] = firstPeakData;
			result[1] = secondPeakData;
		}

	private:
		TGraphErrors fitData;
		TF1 fitModel;
};


class FindDoublePeakWithDynamicPedestal : public FindDoublePeak{
	public:
		// FindDoublePeakWithDynamicPedestal normal constructor with configStruct initilization parameters
		FindDoublePeakWithDynamicPedestal(const configStruct &initParams) : FindDoublePeak(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData)
		{
			fitModel = TF1("fitModel",fitModel,0.0,hitPeriod,5);
			const double timeShift0 = 30.0;
			const double scalingFactor0 = TMath::Max(result[1].peakHeight * bits2scalingFactor, 1000.0);
			const double Q = result[0].peakHeight;
			const double timeshift1 = result[2].peakTime - result[1].peakTime;
			const double scalingFactor1 = TMath::Max(result[2].peakHeight * bits2scalingFactor, 1000.0);

			Double_t initialParameters[5] = {timeShift0, scalingFactor0, Q, timeshift1, scalingFactor1};
			Double_t finalParameters[5];

			adcWaveform2TGraphErrors(adcData, fitData);
			fitModel2Waveform(fitModel, fitData, initialParameters, finalParameters);
			fitParams2ResultantData(result, finalParameters);
		}

	private:
		TGraphErrors fitData;
		TF1 fitModel;


		void fitParams2ResultantData(resultantHitData &result, Double_t *fitParameters)
		{
			resultantPeakData firstPeakData;

			// If there's time maybe make timeshift,scaling factor enumerated like true anomaly and eccentricity
			firstPeakData.peakTime = fitParameters[0];
			firstPeakData.peakHeight = fitParameters[1];

			resultantPeakData secondPeakData;
			secondPeakData.peakTime = fitParameters[3];
			secondPeakData.peakHeight = fitParameters[4];

			result[0] = firstPeakData;
			result[1] = secondPeakData;
			result.pop_back();
		}
};




class FindMultiplePeaks : public FindPeakBase{
	public:

		// FindMultiplePeaks normal constructor with configStruct initilization parameters
		FindMultiplePeaks(const configStruct &initParams) : FindPeakBase(initParams){}

		// Fills result using adc waveform data using by fitting with the convolutionSinglePeakWithDynamicPedestal model
		// NOTE : This function may begin with peak data provided in result which is replaced
		virtual void process(resultantHitData &result, const adcWaveform adcData)
		{
			adcWaveform2TGraphErrors(adcData, fitData);
			findPeaks(fitData,result,initParams);

			const int numPeaks = result.size();

			if (result[0].peakTime == 0.0) // If there is a dynamic pedestal
			{

				// If the only peak is a dynamic pedestal 
				// search for another peak
				if (numPeaks == 1) 
				{
					dynamicPedestalAddPeak(fitData, result);
					FindSinglePeakWithDynamicPedestal singlePeak(initParams);
					singlePeak.process(result, adcData);
				}
				else if (numPeaks == 2)
				{
					FindSinglePeakWithDynamicPedestal singlePeak(initParams);
					singlePeak.process(result, adcData);
				}
				else if (numPeaks == 3)
				{
					FindDoublePeakWithDynamicPedestal doublePeak(initParams);
					doublePeak.process(result, adcData);
				}
			}
			// If there is no dynamic pedestal
			else
			{
				if (numPeaks == 1)
				{
					FindSinglePeak singlePeak(initParams);
					singlePeak.process(result, adcData);
				}
				if (numPeaks == 2)
				{
					FindDoublePeak doublePeak(initParams);
					doublePeak.process(result, adcData);
				}

			}	
		}

	private:
		// Performs explicit peak search on adc waveform data
		void findPeaks(TGraphErrors &gr,resultantHitData &result, const configStruct &initParams, double sigma = 3.0)
		{
  			int ientry = 0; // Start time at 0
  			const double *measurementTimes = gr.GetX();
  			const double *adcValues = gr.GetY();

  			while(ientry < initParams.numSamplesPerHit)
  			{
    			double adcValue = adcValues[ientry];
				double tMax = measurementTimes[ientry];
    			double adcMax = adcValue;
    			double adcPrev = adcValue;

    			int jentry = ientry + 1;
    			bool descending = false;
    			while (jentry < initParams.numSamplesPerHit)
				{
					adcValue = adcValues[jentry];
      				descending |= ((adcPrev-adcValue) > (TMath::Sqrt(2.0)*initParams.adcError*sigma));

      				if (descending && (adcValue-adcPrev > (TMath::Sqrt(2.0)*initParams.adcError*sigma)))
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
		void dynamicPedestalAddPeak(TGraphErrors &gr, resultantHitData &result)
		{	
			// This maybe could be done using linear algebra vectors
			// instead of arrays
			const Double_t *adcValues = gr.GetY();
			const Double_t *measurementTimes = gr.GetX();
			Double_t subtractedValues[initParams.numSamplesPerHit];

			Double_t dynamicPedstalParam[1] = {adcValues[0]};
			Double_t dynamicPedestalX[1];

			for (int i = 0; i < initParams.numSamplesPerHit; ++i)
			{
				dynamicPedestalX[0] = measurementTimes[i];
				subtractedValues[i] = adcValues[i] - dynamicPedestal(dynamicPedestalX, dynamicPedstalParam);
			}

			// New peak is max value of difference between of adc values and dynamic pedestal
			const Float_t newAdcPeak = TMath::MaxElement(initParams.numSamplesPerHit, subtractedValues);
			const Float_t newTPeak = TMath::LocMax(initParams.numSamplesPerHit, subtractedValues);

			resultantPeakData newPeakData(newTPeak, newAdcPeak);
			result.push_back(newPeakData);
		}

		TGraphErrors fitData;
		TF1 fitModel;


};

