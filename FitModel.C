#include "TMath.h"
#include <iostream>
#include "TF1.h"

using namespace std;

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


		float unConvolvedSinglePeak(double *x, double *par)
		{
  		// Initial return value
  		float returnValue = 0.0;
  		// Set x value
  		float xValue = x[0];
  
  		if (xValue > 0.0)
  		{
    	returnValue = (pow(par[0]*xValue/par[1],par[0])
                  /(par[1]*TMath::Gamma(par[0])))
                    *exp(-par[0]*xValue/par[1]);
  		}
  		return returnValue;
		}


		// Note that this is a convolution with a uniform distribution
		//2 Parameters (shaping power set to 1.0)
		//par[0] - Shaping Time
		//par[1] - sigma
		float convolvedSinglePeak(double *x, double *par)
		{
    	float returnValue = 0.0;
		
		if (par[1] == 0.0)
    	{
      		double parameters[2] = {1.0,par[0]}; 
      		returnValue = unConvolvedSinglePeak(x,parameters);
    	}
    	else
    	{
      		const float a = TMath::Max((x[0] + par[1]) / par[0],0.0);
      		// Assuming that shaping time is positive and thus b is negative (if t - sigma is)
      		const float b = TMath::Max((x[0] - par[1]) / par[0],0.0);
      		returnValue =  (-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * par[1]);
    	}
    	return returnValue;
    }

    TF1 *fitFunction;


};


class convolutionSinglePeak : protected FitModelBase{

	public: 

//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is vertical shift 1st peak
//par[3] is sigma 1st peak

		virtual float fitModel(double *x, double *par)
    {

      double convolvedSinglePeakX[1] = {x[0] - par[0]};

      double convolvedSinglePeakParams[2] = {shapingTime,par[3]};

      return (float) par[1] * convolvedSinglePeak(convolvedSinglePeakX,convolvedSinglePeakParams) + par[2];

    }
			
};

// This is a truncating fitting function with a dynamical pedestal
// par[0] is shifted time
// par[1] is scaling factor
// par[2] is Q
// par[3] is sigma

// This should inherit from convolutionSinglePeak not FitModelBase
class convolutionSinglePeakWithDyamicPedestal : FitModelBase{

  public: 
    virtual float fitModel(double *x, double *par)
    {

      double convolutionSinglePeakParams[4] = {par[0],par[1],0.0,par[3]};

      double dynamicPedestalParam[1] = {par[2]};


      // This looks like a awful way of calling the function
      // convolution single peak
      convolutionSinglePeak c;

      return (float) c.convolutionSinglePeak::fitModel(x,convolutionSinglePeakParams) + dynamicPedestal(x,dynamicPedestalParam);
    }
};

