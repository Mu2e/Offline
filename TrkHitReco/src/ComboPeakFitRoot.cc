// fit waveform using root TF1
//
// The member variables of PeakFitRoot need to be protected not private
#include "TrkHitReco/inc/ComboPeakFitRoot.hh"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"

namespace mu2e {

  namespace TrkHitReco {

  ComboPeakFitRoot::ComboPeakFitRoot(const StrawResponse& srep, const fhicl::ParameterSet& pset) : 
      PeakFitRoot(srep, pset) 
  {}



  void ComboPeakFitRoot::process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const {
    // find initial values for the fit
    PeakFit::process(adcData,fit);
    if (_debug>0)std::cout << "PeakFitRoot Initialization charge = " << fit._charge << std::endl;

    // convert waveform to a TGraph
    TGraphErrors fitData;
    adcWaveform2TGraphErrors(adcData,fitData);
    peakResultVector initialGuess;
    findPeaks(fitData, initialGuess, 3.0);
    addEarlyPeak(fitData, initialGuess);

    // set the initial values based on this 'fit'
    Double_t parray[PeakFitParams::nParams];
    fit.fillArray(parray);
    _peakfit.fitModelTF1()->SetParameters(parray);

    int locPrimaryPeak = 0;
    _peakfit.fitModelTF1()->ReleaseParameter(PeakFitParams::width);
    _peakfit.fitModelTF1()->SetParameter(PeakFitParams::width,12.0);

    if (hasEarlyCharge(initialGuess))
    {
        _peakfit.fitModelTF1()->ReleaseParameter(PeakFitParams::earlyCharge);
        _peakfit.fitModelTF1()->SetParameter(PeakFitParams::earlyCharge, initialGuess[0]._peakHeight - (double) _srep.ADCPedestal()); // Are units correct???
        ++locPrimaryPeak;
    }
    // If there should be a floating pedestal
    if(_config.hasOption(FitConfig::floatPedestal))
    {
        _peakfit.fitModelTF1()->ReleaseParameter(PeakFitParams::pedestal);
        _peakfit.fitModelTF1()->SetParameter(PeakFitParams::earlyCharge, _srep.ADCPedestal() ); // Are units correct???
    }
    if (hasLateCharge(initialGuess))
    {
        // Make sure that width is defaulted to being a free parameteri
        _peakfit.fitModelTF1()->FixParameter(PeakFitParams::width, 0.0);
        _peakfit.fitModelTF1()->ReleaseParameter(PeakFitParams::lateCharge);
        _peakfit.fitModelTF1()->SetParameter(PeakFitParams::lateCharge, initialGuess[locPrimaryPeak+1]._peakHeight - _srep.ADCPedestal()); // Are units correct???
        _peakfit.fitModelTF1()->SetParameter(PeakFitParams::lateShift, initialGuess[locPrimaryPeak+1]._peakTime); // Is this shifted correctly
    }
        // debug
        if(_debug>1){
              std::cout << "data = ";
              for (size_t i = 0; i < adcData.size(); ++i){
                  std::cout << fitData.GetY()[i] << ",  ";
                }
              std::cout << std::endl;
              std::cout << "func = ";
                for (size_t i = 0; i < adcData.size(); ++i){
                  std::cout << _peakfit.fitModelTF1()->Eval(fitData.GetX()[i]) << ",  ";
                }
              std::cout << std::endl;
    }

    // invoke the fit
    TFitResultPtr fitresult = fitData.Fit(_peakfit.fitModelTF1(),_fitoptions.c_str());
    // if the fit is a failure, try again, up to the maximum # of iterations
    unsigned ifit=1;
    if(fitresult->Status() ==4 && ifit < _config._maxnit){
      ++ifit;
      fitData.Fit(_peakfit.fitModelTF1(),_fitoptions.c_str());
    }
    // copy fit result back to the the PeakFitParams
    fit = PeakFitParams(_peakfit.fitModelTF1()->GetParameters(),
      fitresult->Chi2(),
      fitresult->Ndf(),
      fitresult->Status());
  }

    bool ComboPeakFitRoot::hasEarlyCharge(const peakResultVector &initialGuess) const
    {
        // Is this stable??
        return initialGuess[0]._peakTime < (_srep.nADCPreSamples()*_srep.adcPeriod());
    }

    bool ComboPeakFitRoot::hasLateCharge(const peakResultVector &initialGuess) const
    {
        // If there are more than 2 peaks or there is more than 1 peak but no early extra peak return true
        int nPeaks = 0; // number of peaks excluding those which appear in the presample data
        for (auto peak : initialGuess)
        {
            if (peak._peakTime >= (_srep.nADCPreSamples()*_srep.adcPeriod())) ++nPeaks;
        }
        return nPeaks>1;
    }


    // Currently the location of late charge is computed using the locLateCharge variable above
    // This is not good design and should be replaced by a function like this
/**    int ComboPeakFitRoot::locLateCharge(const peakResultVector initialGuess) const
    {
        int numEarlyCharge = 0;
        while(hasEarlyCharge(initialGuess))
        {
            initialGuess.erase(initialGuess.begin());
            ++numEarlyCharge;
        }
        return numEarlyCharge+1;
    }**/

    // TODO : adcErrors needs to be gotten from strawele


    //TODO: Convert _initParams._numSamplesPerHit to corresponding element of strawele

    void ComboPeakFitRoot::addEarlyPeak(const TGraphErrors &gr, peakResultVector &initialGuess) const
    {
      //This maybe could be done using linear algebra vectors
      //instead of arrays
      const Double_t *adcValues = gr.GetY();
      const Double_t *measurementTimes = gr.GetX();
      const int numSamplesPerHit = gr.GetN();
      Double_t subtractedValues[numSamplesPerHit];

      const double earlyPeakCharge = adcValues[0];

      for (int i = 0; i < numSamplesPerHit; ++i)
      {
         PeakFitFunction func(_srep);
         func.init(_config);
         // FIXME : THIS NEEDS TO GET THE FUNCTION EARLYPEAK FROM PEAK FIT FUNCTION
         subtractedValues[i] = adcValues[i] - func.earlyPeak(measurementTimes[i], earlyPeakCharge);
      }

      //New peak is max value of difference between of adc values and dynamic pedestal
      const Float_t newAdcPeak = TMath::MaxElement(numSamplesPerHit, subtractedValues);
      const Float_t newTPeak = TMath::LocMax(numSamplesPerHit, subtractedValues);

      peakResult newPeakData(newTPeak, newAdcPeak);
      initialGuess.push_back(newPeakData);
    }

    // Performs explicit peak search on adc waveform data
    void ComboPeakFitRoot::findPeaks(const TGraphErrors &gr, peakResultVector &initialGuess, const double sigma) const
    {
      int ientry = 0; // Start time at 0
      const double *measurementTimes = gr.GetX();
      const double *adcValues = gr.GetY();
      const int numSamplesPerHit = gr.GetN();

      while(ientry < numSamplesPerHit)
      {
        double adcValue = adcValues[ientry];
        double tMax = measurementTimes[ientry];
        double adcMax = adcValue;
        double adcPrev = adcValue;

        int jentry = ientry + 1;
        bool descending = false;
        while (jentry < numSamplesPerHit)
        {
          adcValue = adcValues[jentry];
          descending |= ((adcPrev-adcValue) > (TMath::Sqrt2()*_srep.analogNoise(StrawElectronics::adc)/_srep.adcLSB()*sigma));

          if (descending && (adcValue-adcPrev > (TMath::Sqrt2()*_srep.analogNoise(StrawElectronics::adc)/_srep.adcLSB()*sigma)))
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
        peakResult peakData(tMax, adcMax - _srep.ADCPedestal());
        initialGuess.push_back(peakData);
        ++ientry;
      }
    }

  } // TrkHitReco namespace

}// mu2e namespace
