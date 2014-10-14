#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <sstream>
#include <vector>

// Par0 - Q
// Par1 - shapingTime
// The dynamic pedestal is of the form Q e^(-t / tau). 
//Note that normal it is Q / tau rather than Q 
float dynamicPedestal(double *x, double *par)
{
    return par[0] * exp(-x[0] / par[1]);

}

// Current function used to calculate current
float currentFunction(double *x, double *par)
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

//2 Parameters. Note that shaping power is fixed to 1
//par[0]- Shaping time 
//par[1] - sigma

float currentFunction2(double *x, double *par)
{
  // Initial return value
  Float_t returnValue = 0.0;
  // Set x value
  float y = x[0];

  float r = par[0];
  float a = par[1];

  Float_t zero = 0.0;

  returnValue =  -1 / (2 * pow(r,3)) * exp(-y/r + a*a / (2 * r * r)) * 
        (a*a - r * y ) * (1 - TMath::Erf((a*a - r * y) / (sqrt(2) * a * r))) + 
        a / (sqrt(2 * TMath::Pi()) * r*r) * exp(-y*y/(2 * a * a));  


  return TMath::Max(returnValue, zero);
}



//Fitting function for current Function2
//par[0] is shifted time 1st peak
//par[1] is scalingfactor 1st peak
//par[2] is vertical shift 1st peak
//par[3] is sigma 1st peak
//par[4] is the level of truncation (set to either 1023 or 1023 - 64)
//par[5] is the shaping time
//shaping time is set to 100.0


// Note that this function truncates above 1023-64 bits
float fittingFunction4(double *x, double *par)
{
    double currentX[1] = {x[0] - par[0]};
    double currentParameters[2] = {par[5],par[3]};

    double truncatingValue = par[4];

    return TMath::Min((par[1] * currentFunction2(currentX, currentParameters) + par[2]),truncatingValue);
}


// Sum of two convolved fitting functions (fittingFunction4)
//par[0] is shifted time
//par[1] is scalingfactor
//par[2] is vertical shift
//par[3] is sigma
float fittingFunction5(double *x, double *par)
{
    double firstParam[4] = {par[0], par[1], par[2], par[3]};
    double secondParam[4] = {par[0] + par[4], par[5], par[6], par[7]};

    return fittingFunction4(x,firstParam)
          +fittingFunction4(x, secondParam);
}


// Note that this function DOES NOT truncate
float fittingFunction6(double *x, double *par)
{
    double currentX[1] = {x[0] - par[0]};
    double currentParameters[2] = {100.0,par[3]};

    return (par[1] * currentFunction2(currentX, currentParameters) + par[2]);
}



// This is a truncating fitting function with a dynamical pedestal
// par[0] is shifted time
// par[1] is scaling factor
// par[2] is Q
// par[3] is sigma
// par[4] is shapin time
// Note that this function depends on fittingFunction4
float fittingFunction7(double *x, double *par)
{
    // vertical shift is set to 0.0 so that dynamic pedestal can be used
    double fittingFunction4Parameters[6] = {par[0],par[1],0.0,par[3], 1023.0 - 64.0,par[4]};
    
    // Shaping time is set to 100.0
    double dynamicPedestalParameters[2] = {par[2],par[4]};

    return fittingFunction4(x,fittingFunction4Parameters) + dynamicPedestal(x,dynamicPedestalParameters);

}




// Fitting funtion used for fitting
// par[0] is shifted time
// par[1] is scaling factor
// par[2] is constant pedestal
// par[3] is shaping time
float fittingFunction1(double *x, double *par)
{
    // These values must be doubles here

    //Shaping time is now free parameter
    //double shapingTime = 100.0;
    double shapingPower = 1.0;

    double currentX[1] = {x[0] - par[0]};
    double currentParameters[2] = {shapingPower,par[3]};

    return (par[1] * currentFunction(currentX, currentParameters)) + par[2];
}


// Sum of two fittingFunctions
float fittingFunction2(double *x, double *par)
{

    double firstParam[3] = {par[0], par[1], par[2]};
    double secondParam[3] = {par[3], par[4], par[5]};

    //double firstX[1] = {x[0]};
   // double secondX[1] = {x[1]};

    return fittingFunction1(x,firstParam)
          +fittingFunction1(x, secondParam);
}

float fittingFunction3(double *x, double *par)
{
    // Par0 - shift in X 1st peak
    // Par1 - scalingFactor 1st peak
    // Par2 - vertical shift 1st peak
    // Par3 - shift in 2nd peak minus shift in 1st peak
    // Par4 - scaling factor 2nd peak
    // Par5 - vertical shift 2nd peak
    // Par6 - shaping time
    double firstParam[4] = {par[0], par[1], par[2], par[6]};
    double secondParam[4] = {par[3] + par[0], par[4], par[5], par[6]};

    //double firstX[1] = {x[0]};
   // double secondX[1] = {x[1]};

    return fittingFunction1(x, firstParam)
          +fittingFunction1(x, secondParam);

}

// This function fits using a double peak and a dynamic pedestal
// Depends on fitting function. Note that the fits are not convoluted? convolved?
// Par0 - shift in X 1st peak
// Par1 - scalingFactor 1st peak
// Par2 - Q
// Par3 - shift in 2nd peak minus shift in 1st peak
// Par4 - scaling factor 2nd peak
// Par5 - shaping time
float fittingFunction8(double *x, double *par)
{
    // Note that pedestals are set to 0.0
    double doublePeakParam[7] = {par[0],par[1],0.0,par[3],par[4],0.0,par[5]};

    // Shaping time is set to 100.0
    double dynamicPedestalParameters[2] = {par[2],par[5]};

    // 
    float truncationValue = 1023.0 - 64.0;

    return TMath::Min(fittingFunction3(x,doublePeakParam) + dynamicPedestal(x,dynamicPedestalParameters), truncationValue);
}

// Just like fitting function 9 but fixes sigma at 7 instead of 0
/**float fittingFunction9(double *x, double *par)
{
    // Note that pedestals are set to 0.0
    double doublePeakParam[7] = {par[0],par[1],0.0,par[3],par[4],0.0,par[5]};

    // Shaping time is set to 100.0
    double dynamicPedestalParameters[2] = {par[2],par[5]};

    // 
    float truncationValue = 1023.0 - 64.0;

    return TMath::Min(fittingFunction3(x,doublePeakParam) + dynamicPedestal(x,dynamicPedestalParameters), truncationValue);
} **/



void convert2StringInt(TString &string, double doubleNum)
{
  int num = (int) doubleNum;
  ostringstream convert;
  convert << num; 
  string = convert.str();
}

void convert2StringInt(TString &string, int intNum)
{
  ostringstream convert;
  convert << intNum; 
  string = convert.str();
}

void convert2StringDouble(TString &string, double doubleNum)
{
  ostringstream convert;
  convert << doubleNum; 
  string = convert.str();
}

TGraph* computeRejectionGraph(TH1F *electronHist, TH1F *protonHist, const int numberOfBins)
{
  Double_t truncX[numberOfBins], truncY[numberOfBins];

  int protonSum = 0;
  int electronSum = 0;
  for (int i = 1; i <= numberOfBins; ++i)
  {
    // For some reason bin number starts with 1 for TH1F
    protonSum += protonHist->GetBinContent(i);
    electronSum += electronHist->GetBinContent(i);
    // acceptance rate of electrons
    truncX[i - 1] = electronSum / (double) electronHist->GetEntries();
    // 1 - rejection rate
    truncY[i - 1] = 1 - (protonSum / (double) protonHist->GetEntries());
  }

  TGraph * rejectionGraph = new TGraph(numberOfBins,truncX,truncY);
  return rejectionGraph;

}



void findPeaks(TGraphErrors *gr,vector<Float_t>& tPeak, vector<Float_t>& adcPeak, double sigma = 3.0)
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

