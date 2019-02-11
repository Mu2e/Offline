#include "TF1.h"
#include "TrkHitReco/inc/PeakFitFunction.hh"
#include <algorithm>
#include <iostream>

namespace mu2e {

  namespace TrkHitReco {

    PeakFitFunction::PeakFitFunction(const StrawResponse& srep) : 
      _srep(srep), _fitConfig(), _tf1(0)
    {}
    PeakFitFunction::~PeakFitFunction() { delete _tf1; }
    
    
    //initialize the fit object after complete config, 
    //NB that the default values of this TF1 must be initialized before fitting!!!!!
    void PeakFitFunction::init(const FitConfig& fitConfig)      
    {
        _fitConfig=fitConfig;
        createTF1();
    }
    
    // the actual fit function, taking the explicit parameters as input
    Double_t PeakFitFunction::fitModel(Double_t time, PeakFitParams const& params) const {
      // start with the pedestal value
      Float_t returnValue(params._pedestal);
      // add the main peak
      returnValue += params._charge*convolvedSinglePeak(time-params._time,params._width);
      // if requested, add the early and late peaks
      if(_fitConfig.hasOption(FitConfig::earlyPeak)){
	// note this is WRT absolute time, not shift time (???)
	returnValue += earlyPeak(time,params._earlyCharge);
      }
      if(_fitConfig.hasOption(FitConfig::latePeak)){
	// note this is WRT absolute time, not shift time (???)
	returnValue += params._lateCharge*convolvedSinglePeak(time-params._time-params._lateShift,params._width);
      }
      // optionally truncate the final result
      if(_fitConfig.hasOption(FitConfig::truncateADC)){
	truncateResponse(returnValue);
      }

      return returnValue;
    }

    void PeakFitFunction::resetTF1(TF1& tf1) { (*_tf1) = tf1; }

    void PeakFitFunction::resetTF1(PeakFitParamsLimits const& params) {
      
      Double_t array[PeakFitParams::nParams];
      params.fillArray(array);
      _tf1->SetParameters(array);

      for (unsigned ipar=0;ipar< PeakFitParams::nParams;++ipar)
      {
	 PeakFitParams::paramIndex ip = static_cast<PeakFitParams::paramIndex>(ipar);
	 if(params.isFree(ip)){
	   _tf1->ReleaseParameter(ipar);
	 } else {
	   _tf1->FixParameter(ipar,_tf1->GetParameter(ipar));
         }
	 
	 _tf1->SetParLimits(ipar,params._parmin[ipar], params._parmax[ipar]);
      }
    }


    Double_t PeakFitFunction::operator()(Double_t* x, Double_t* p) {
      return fitModelRoot(x,p); 
    }

    // the root version of same.  This calls down to the above
    Double_t PeakFitFunction::fitModelRoot(Double_t* x, Double_t* p) const {
      Double_t time = x[0];
      PeakFitParams params(p);
      Float_t result = fitModel(time,params);
//      std::cout << "PeakFitFunction time = " << time 
//	<< " peak time = " << params._time 
//	<< " peak charge = " << params._charge 
//	<< " model = " << result << std::endl;
      return result;
    }


    // Method for creating a TF1 using the fit function.  This OVERWRITES the state of the TF1
    void PeakFitFunction::createTF1() {
      // bind the location of the root-syntax C-style function
      FitFunctionRoot rootfun = std::bind(&PeakFitFunction::fitModelRoot,this,std::placeholders::_1,std::placeholders::_2);
      // create the TF1
      // NB: the min and max times should come from a global config object, FIXME!!!!
//      _tf1 = new TF1("TrkHitReco::PeakFitFunction",rootfun,0.0,1695.0,PeakFitParams::nParams);
       _tf1 = new TF1("TrkHitReco::PeakFitFunction",*this,0.0,1695.0,PeakFitParams::nParams);
      // Set the parameter names
      for(size_t iparam=0;iparam< PeakFitParams::nParams; ++ iparam)
	_tf1->SetParName(iparam,PeakFitParams::parameterName((PeakFitParams::paramIndex) iparam).c_str());
     // decide which parameters are fixed and which are free, and initialize the dormant parameters
      //
      //  No early peak: fix the early charge to 0
      if(!_fitConfig.hasOption(FitConfig::earlyPeak))
	_tf1->FixParameter(PeakFitParams::earlyCharge,0.0);
      //no floating pedestal: fix the value to what strawelectronics says
      if(!_fitConfig.hasOption(FitConfig::floatPedestal))
	_tf1->FixParameter(PeakFitParams::pedestal,_srep.ADCPedestal());
      // no floating width: fix the value to what strawelectronics says
      if(!_fitConfig.hasOption(FitConfig::floatWidth))
	_tf1->FixParameter(PeakFitParams::width,_srep.fallTime(StrawElectronics::adc));
      // no late peak: fix the charge and shift and disable those parameters
      if(!_fitConfig.hasOption(FitConfig::latePeak)){
	_tf1->FixParameter(PeakFitParams::lateShift,0.0);
	_tf1->FixParameter(PeakFitParams::lateCharge,0.0);
      }
      // time and charge (of the main peak) are ALWAYS left as free parameters
      // set parameter errors and limits
      _tf1->SetParError(PeakFitParams::charge,0.01);
      _tf1->SetParError(PeakFitParams::time,5.0);
      _tf1->SetParError(PeakFitParams::earlyCharge,0.1);
      _tf1->SetParLimits(PeakFitParams::charge,0.0,10.0);
      _tf1->SetParLimits(PeakFitParams::time,0.0,80.0);
      _tf1->SetParLimits(PeakFitParams::earlyCharge,0.0,1.0);
      _tf1->SetParLimits(PeakFitParams::lateShift,20.0,70.0);
      _tf1->SetParLimits(PeakFitParams::lateCharge,0.0,10.0);
      // limit the width to be > 0
      _tf1->SetParLimits(PeakFitParams::width,0.0,30.0);
      // limit the pedestal to +- 5 sigma noise
      double pednoise =_srep.analogNoise(StrawElectronics::adc)/_srep.adcLSB();
      double pedmin = std::max(0.0,_srep.ADCPedestal()-5.0*pednoise);
      double pedmax = _srep.ADCPedestal()+5.0*pednoise;
      _tf1->SetParLimits(PeakFitParams::pedestal,pedmin,pedmax);
    }

    // model charge draining from an earlier hit.  Only the exponential portion should remain
    // need to check/fix the normalization, FIXME!!!
    Float_t PeakFitFunction::earlyPeak(const Double_t time, const Double_t charge) const
    {
      return charge * exp(-time / _srep.fallTime(StrawElectronics::adc));
    }

    // Normalized CR-RC network response to a current delta function
    Float_t PeakFitFunction::unConvolvedSinglePeak(const Double_t time) const
    {
       Float_t returnValue = 0.0;
       if (time > 0.0)
       {
	 const double pC_per_uA_ns{1000}; // unit conversion from pC/ns to microAmp
	 static const double norm = _srep.currentToVoltage(StrawElectronics::adc)*
	 pow(_srep.fallTime(StrawElectronics::adc),-2)/pC_per_uA_ns;
	 returnValue = time*norm*exp(-time/_srep.fallTime(StrawElectronics::adc));
       }
       return returnValue;
    }

    
    // Same as above, but convoluted with a uniform distribution
    // need to fix the normalization FIXME!!!
    Float_t PeakFitFunction::convolvedSinglePeak(const Double_t time, const Double_t sigma) const
    {
	const double pC_per_uA_ns{1000}; // unit conversion from pC/ns to microAmp
        static const double norm = _srep.currentToVoltage(StrawElectronics::adc)/pC_per_uA_ns;
        Float_t returnValue = 0.0;

        if (sigma <= 0.0)
        {
	  returnValue = unConvolvedSinglePeak(time);
        }
        else
        {
	  const Float_t a = std::max((time + sigma) / _srep.fallTime(StrawElectronics::adc),0.0);
	  // Assuming that shaping time is pggositive and thus b is negative (if t - sigma is)
	  const Float_t b = std::max((time - sigma) / _srep.fallTime(StrawElectronics::adc),0.0);
	  returnValue =  norm*(-exp(-a)*(1+a) + exp(-b)*(1+b)) / (2.0 * sigma);
	  // this value doesn't have the correct absolute normalization, FIXME!!!!
        }
        return returnValue;
    }
    
    
    
    //  model the truncation.  We first have to translate from ADC units to voltage, then back (!)
    void PeakFitFunction::truncateResponse(Float_t& rval) const 
    {
      double mv = _srep.adcLSB()*(rval-_srep.ADCPedestal());
      double saturatedmv = _srep.saturatedResponse(mv);
      double saturatedadc = saturatedmv/_srep.adcLSB() + _srep.ADCPedestal();
      static const double maxadc(_srep.maxADC());
      rval = std::max(std::min(saturatedadc,maxadc),(double)0.0);
    }
  }
}
