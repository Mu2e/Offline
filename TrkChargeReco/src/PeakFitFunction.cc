#include "TF1.h"
#include "TrkChargeReco/inc/PeakFitFunction.hh"
#include <algorithm>


namespace mu2e {

  namespace TrkChargeReco {
    // constructor
    PeakFitFunction::PeakFitFunction(StrawElectronics const& strawele, FitConfig const& fitConfig) : _strawele(strawele),
    _fitConfig(fitConfig), _tf1(0) {
      // create the TF1 for root fitting.  NB that the default values of this TF1 must be initialized before fitting!!!!!
      createTF1();
    }
    PeakFitFunction::~PeakFitFunction() { delete _tf1; }
    // the actual fit function, taking the explicit parameters as input
    Float_t 
      PeakFitFunction::fitModel(Double_t time, PeakFitParams const& params) const {
	// start with the pedestal value
	Float_t returnValue(params._pedestal);
	// add the main peak
	returnValue += params._scale*convolvedSinglePeak(time-params._time,params._width);
	// if requested, add the early and late peaks
	if(_fitConfig.hasOption(FitConfig::earlyPeak)){
	  // note this is WRT absolute time, not shift time (???)
	  returnValue += earlyPeak(time,params._earlyCharge);
	}
	if(_fitConfig.hasOption(FitConfig::latePeak)){
	  // note this is WRT absolute time, not shift time (???)
	  returnValue += params._lateScale*convolvedSinglePeak(time-params._time-params._lateShift,params._width);
	}
	// optionally truncate the final result
	if(_fitConfig.hasOption(FitConfig::truncateADC)){
	  truncateResponse(returnValue);
	}

	return returnValue;
      }

    // the root version of same.  This calls down to the above
    Float_t PeakFitFunction::fitModelRoot(Double_t* x, Double_t* p) const {
      // copy the input parameters
      Double_t time = x[0];
      PeakFitParams params(p);
      return fitModel(time,params);
    }
    // Method for creating a TF1 using the fit function.  This OVERWRITES the state of the TF1
    void PeakFitFunction::createTF1() {
      // bind the location of the root-syntax C-style function
      FitFunctionRoot rootfun = std::bind(&PeakFitFunction::fitModelRoot,this,std::placeholders::_1,std::placeholders::_2);
      // create the TF1
      // NB: the min and max times should come from a global config object, FIXME!!!!
      _tf1 = new TF1("TrkChargeReco::PeakFitFunction",rootfun,0.0,1695.0,PeakFitParams::nParams);
      // decide which parameters are fixed and which are free, and initialize the dormant parameters
      //
      //  No early peak: fix the early charge to 0
      if(!_fitConfig.hasOption(FitConfig::earlyPeak))
	_tf1->FixParameter(PeakFitParams::earlyCharge,0.0);
      // no floating pedestal: fix the value to what strawelectronics says
      if(!_fitConfig.hasOption(FitConfig::floatPedestal))
	_tf1->FixParameter(PeakFitParams::pedestal,_strawele.ADCPedestal());
      // no floating width: fix the value to what strawelectronics says
      if(!_fitConfig.hasOption(FitConfig::floatWidth))
	_tf1->FixParameter(PeakFitParams::width,_strawele.fallTime(StrawElectronics::adc));
      // no late peak: fix the scale and shift and disable those parameters
      if(!_fitConfig.hasOption(FitConfig::latePeak)){
	_tf1->FixParameter(PeakFitParams::lateShift,0.0);
	_tf1->FixParameter(PeakFitParams::lateScale,0.0);
      }
      // time and scale (of the main peak) are ALWAYS left as free parameters
    }

    // model charge draining from an earlier hit.  Only the exponential portion should remain
    Float_t PeakFitFunction::earlyPeak(const Double_t time, const Double_t charge) const
    {
      return charge * exp(-time / _strawele.fallTime(StrawElectronics::adc));
    }

    // Normalized CR-RC network response to a current delta function
    Float_t PeakFitFunction::unConvolvedSinglePeak(const Double_t time) const
    {
      // Initial return value
      Float_t returnValue = 0.0;

      if (time > 0.0)
      {
	static double norm = pow(_strawele.fallTime(StrawElectronics::adc),-2);
	returnValue = time*norm*exp(-time/_strawele.fallTime(StrawElectronics::adc));
      }
      return returnValue;
    }

    // Same as above, but convoluted with a uniform distribution
    Float_t PeakFitFunction::convolvedSinglePeak(const Double_t time, const Double_t sigma) const
    {
      Float_t returnValue = 0.0;

      if (sigma <= 0.0)
      {
	returnValue = unConvolvedSinglePeak(time);
      }
      else
      {
	const Float_t a = std::max((time + sigma) / _strawele.fallTime(StrawElectronics::adc),0.0);
	// Assuming that shaping time is pggositive and thus b is negative (if t - sigma is)
	const Float_t b = std::max((time - sigma) / _strawele.fallTime(StrawElectronics::adc),0.0);
	returnValue =  (-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * sigma);
      }
      return returnValue;
    }
    //  model the truncation.  We first have to translate from ADC units to voltage, then back (!)
    void PeakFitFunction::truncateResponse(Float_t& rval) const {
      double mv = _strawele.adcLSB()*(rval-_strawele.ADCPedestal());
      double saturatedmv = _strawele.saturatedResponse(mv);
      double saturatedadc = saturatedmv/_strawele.adcLSB() + _strawele.ADCPedestal();
      static const double maxadc(_strawele.maxADC());
      rval = std::max(std::min(saturatedadc,maxadc),(double)0.0);
    }
  }
}
