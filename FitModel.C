#include "TMath.h"
#include <iostream>
#include "TF1.h"

using namespace std;

// ADD FIXED TRUNCATION
 
class FitModelBase{
	public:
		virtual float fitModel(double *x, double *par) = 0;
		virtual ~FitModelBase(){}

    double shapingTime = 25.0;
    double truncationLevel = 1023.0 - 64.0;

    // Par0 - Q
// TMath::Min((par[1 - shapingTime
// The dynamic pedestal is of the form Q e^(-t / tau). 
//Note that normalized it is Q / tau rather than Q 
    float dynamicPedestal(double *x, double *par)
    {
    return par[0] * exp(-x[0] / par[1]);
    }


	protected:

    float fixedTruncation(float currentFunctionValue)
    {
        return TMath::Min(currentFunctionValue,(float)(1023.0-64.0));
    }


    // Shaping power set to 1
    // MAYBE GET RID OF PAR
		float unConvolvedSinglePeak(double *x, double *par)
		{
  		// Initial return value
  		float returnValue = 0.0;
  		// Set x value
  		float xValue = x[0];
  
  		if (xValue > 0.0)
  		{
    	returnValue = xValue*pow(shapingTime,-2)
                    *exp(-xValue/shapingTime);
  		}
  		return returnValue;
		}


		// Note that this is a convolution with a uniform distribution
		//2 Parameters (shaping power set to 1.0)
		//par[0] - sigma
		float convolvedSinglePeak(double *x, double *par)
		{
    	float returnValue = 0.0;
		
		if (par[0] == 0.0)
    	{
      		double parameters[2] = {1.0,shapingTime}; 
      		returnValue = unConvolvedSinglePeak(x,parameters);
    	}
    	else
    	{
      		const float a = TMath::Max((x[0] + par[0]) / shapingTime,0.0);
      		// Assuming that shaping time is positive and thus b is negative (if t - sigma is)
      		const float b = TMath::Max((x[0] - par[0]) / shapingTime,0.0);
      		returnValue =  (-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * par[0]);
    	}
    	return returnValue;
    }

    TF1 *fitFunction;


};


class convolutionSinglePeak : protected FitModelBase{

	public: 

  public:
//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is sigma 1st peak
		virtual float fitModel(double *x, double *par)
    {

      double convolvedSinglePeakX[1] = {x[0] - par[0]};

      double convolvedSinglePeakParams[2] = {par[2]};

      return (float) par[1] * convolvedSinglePeak(convolvedSinglePeakX,convolvedSinglePeakParams);

    }
			
};

class convolutionSinglePeakWithConstantPedestal : protected convolutionSinglePeak{

public:
//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is vertical shift 1st peak
//par[3] is sigma 1st peak

    virtual float fitModel(double *x, double *par)
    {
      double convolutionSinglePeakParams[3] = {par[0],par[1],par[3]};

      return (float) convolutionSinglePeak::fitModel(x,convolutionSinglePeakParams) + par[2];
    }

};

// This is a truncating fitting function with a dynamical pedestal
// par[0] is shifted time
// par[1] is scaling factor
// par[2] is Q
// par[3] is sigma

// This should inherit from convolutionSinglePeak not FitModelBase
class convolutionSinglePeakWithDynamicPedestal : protected convolutionSinglePeak{

  public: 
    virtual float fitModel(double *x, double *par)
    {

      double convolutionSinglePeakParams[4] = {par[0],par[1],0.0,par[3]};

      double dynamicPedestalParam[1] = {par[2]};

      return (float) convolutionSinglePeak::fitModel(x,convolutionSinglePeakParams) 
                    + FitModelBase::dynamicPedestal(x,dynamicPedestalParam);
    }
};

// PUT THE DOUBLE PEAK FUNCTIONS HERE
class doublePeak : public FitModelBase{
  public:
    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - shift in 2nd peak minus shift in 1st peak
    // Par3 - scaling factor 2nd peak
    virtual float fitModel(double *x, double *par)
    {
      // Why can't these be constant arrays???
      convolutionSinglePeak c;
      double firstPeakPar[3] = {par[0], par[1], 0.0};
      double secondPeakPar[3] = {par[2] + par[0], par[3], 0.0};

      return c.fitModel(x, firstPeakPar) 
           + c.fitModel(x, secondPeakPar);
    }

};

class doublePeakWithConstantPedestal : public doublePeak{
public:

  // Par0 - shift in X 1st peak
  // Par1 - scalingFactor 1st peak
  // Par2 - vertical shift
  // Par3 - shift in 2nd peak minus shift in 1st peak
  // Par4 - scaling factor 2nd peak
  virtual float fitModel(double *x, double *par)
  {
    double doublePeakParams[4] = {par[0],par[1],par[3],par[4]};
    return doublePeak::fitModel(x, doublePeakParams) + par[2];
  }

};

// Par0 - shift in X 1st peak
// Par1 - scalingFactor 1st peak
// Par2 - Q
// Par3 - shift in 2nd peak minus shift in 1st peak
// Par4 - scaling factor 2nd peak
class doublePeakWithDynamicPedestal : public doublePeak{
public:
  virtual float fitModel(double *x, double *par)
  {
    double doublePeakParams[4] = {par[0],par[1],par[3],par[4]};
    double dynamicPedestalParam[1] = {par[2]};

    return doublePeak::fitModel(x,doublePeakParams) 
        + FitModelBase::dynamicPedestal(x,dynamicPedestalParam);

  }

};
