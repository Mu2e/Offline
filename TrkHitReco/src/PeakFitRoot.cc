// fit waveform using root TF1
#include "TrkHitReco/inc/PeakFitRoot.hh"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include <iostream>

namespace mu2e {
  
  namespace TrkHitReco {

    PeakFitRoot::PeakFitRoot(const StrawResponse& srep, const fhicl::ParameterSet& pset) : 
        PeakFit(srep,pset),
        _truncateADC(pset.get<bool>(      "TruncateADC",true)), 
        _floatPedestal(pset.get<bool>(    "FloatPedestal",true)), 
        _floatWidth(pset.get<bool>(       "FloatWidth",true)), 
        _earlyPeak(pset.get<bool>(        "EarlyPeak",false)),
        _latePeak(pset.get<bool>(         "LatePeak",false)),
        _fitoptions(pset.get<std::string>("fitOption","QNSEX0B")),
        _maxFitIter(pset.get<unsigned>(   "MaxFitIterations",1)),
        _debug(pset.get<int>(             "debugLevel",0)),
        _peakfit(srep),
        _config(_maxFitIter, _debug)
    {   
        if (_floatWidth)    _config.setOption(TrkHitReco::FitConfig::floatWidth);
        if (_floatPedestal) _config.setOption(TrkHitReco::FitConfig::floatPedestal);
        if (_truncateADC)   _config.setOption(TrkHitReco::FitConfig::truncateADC);
        if (_earlyPeak)     _config.setOption(TrkHitReco::FitConfig::earlyPeak);
        if (_latePeak)      _config.setOption(TrkHitReco::FitConfig::latePeak);
        _peakfit.init(_config);
    }


    void PeakFitRoot::process(TrkTypes::ADCWaveform const& adcData, PeakFitParams & fit) const 
    {

       // find initial values for the fit and set the initial values based on this 'fit'
       PeakFit::process(adcData,fit);
       if (_debug>0)std::cout << "PeakFitRoot Initialization charge = " << fit._charge << std::endl;
       Double_t parray[PeakFitParams::nParams];
       fit.fillArray(parray);
       _peakfit.fitModelTF1()->SetParameters(parray);

       // temporarily set width negative to force an unconvolved fit FIXME!!
       _peakfit.fitModelTF1()->SetParameter(PeakFitParams::width,-1.0);

       TGraphErrors fitData;
       adcWaveform2TGraphErrors(adcData,fitData);

       if (_debug>1)
       {
         std::cout << "data = ";
         for (size_t i = 0; i < adcData.size(); ++i)
	   std::cout << fitData.GetY()[i] << ",  ";
         std::cout << std::endl;
         std::cout << "func = ";  
         for (size_t i = 0; i < adcData.size(); ++i)
	   std::cout << _peakfit.fitModelTF1()->Eval(fitData.GetX()[i]) << ",  ";
         std::cout << std::endl;
       }

       // invoke the fit
       TFitResultPtr fitresult = fitData.Fit(_peakfit.fitModelTF1(),_fitoptions.c_str()); 

       unsigned ifit=1;
       if (fitresult->Status()==4 && ifit <_maxFitIter)
       {
         ++ifit;
         fitData.Fit(_peakfit.fitModelTF1(),_fitoptions.c_str());
       }

       fit = PeakFitParams(_peakfit.fitModelTF1()->GetParameters(),
       fitresult->Chi2(),
       fitresult->Ndf(),
       fitresult->Status());

       for(int ipar=0;ipar < _peakfit.fitModelTF1()->GetNpar();++ipar)
       {
          Double_t parmin, parmax;
          _peakfit.fitModelTF1()->GetParLimits(ipar,parmin,parmax);
          if(parmin == parmax)
	    fit.fixParam((TrkHitReco::PeakFitParams::paramIndex)ipar);
          else
	    fit.freeParam((TrkHitReco::PeakFitParams::paramIndex)ipar);
       }
    }

    void PeakFitRoot::adcWaveform2TGraphErrors(TrkTypes::ADCWaveform const& adcData, TGraphErrors &fitData) const
    {
        Double_t adcDataTemp[adcData.size()];
        Double_t measurementTimes[adcData.size()];
        Double_t measurementTimesErrors[adcData.size()];
        Double_t adcDataErrors[adcData.size()];

        for (size_t i = 0; i < adcData.size(); ++i)
        {
	    adcDataTemp[i] = (Double_t) adcData[i];
	    measurementTimes[i] = (Double_t) i * _srep.adcPeriod(); // should deal with global time offset FIXME!!
	    measurementTimesErrors[i] = 0.0; //_srep.adcPeriod();
	    adcDataErrors[i] = 1.0*_srep.analogNoise(StrawElectronics::adc)/_srep.adcLSB(); // should be able to scale the error FIXME!!
        }
        fitData = TGraphErrors(adcData.size(),measurementTimes,adcDataTemp,measurementTimesErrors,adcDataErrors);
     }


  } 
}


