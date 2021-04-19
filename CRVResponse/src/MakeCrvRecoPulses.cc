#include "CRVResponse/inc/MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>

namespace
{
  double Gumbel(double* xs, double* par)
  {
    double const x = xs[0];
    return par[0]*(TMath::Exp(-(x-par[1])/par[2]-TMath::Exp(-(x-par[1])/par[2])));
  }
  double Gumbel2(double* xs, double* par)
  {
    double const x = xs[0];
    double y1=par[0]*(TMath::Exp(-(x-par[1])/par[2]-TMath::Exp(-(x-par[1])/par[2])));
    double y2=par[3]*(TMath::Exp(-(x-par[4])/par[5]-TMath::Exp(-(x-par[4])/par[5])));
    return y1+y2;
  }
}

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(float minADCdifference, float defaultBeta, float minBeta, float maxBeta,
                                     float maxTimeDifference, float minPulseHeightRatio, float maxPulseHeightRatio,
                                     float LEtimeFactor, bool allowDoubleGumbel, float doubleGumbelThreshold) :
                                     _f1("peakfitter",Gumbel,0,0,3), _f2("peakfitter",Gumbel2,0,0,6),
                                     _minADCdifference(minADCdifference), 
                                     _defaultBeta(defaultBeta), _minBeta(minBeta), _maxBeta(maxBeta), 
                                     _maxTimeDifference(maxTimeDifference),
                                     _minPulseHeightRatio(minPulseHeightRatio),
                                     _maxPulseHeightRatio(maxPulseHeightRatio),
                                     _LEtimeFactor(LEtimeFactor),
                                     _allowDoubleGumbel(allowDoubleGumbel), _doubleGumbelThreshold(doubleGumbelThreshold)
{}

void MakeCrvRecoPulses::FillGraphAndFindPeaks(const std::vector<unsigned int> &waveform, unsigned int startTDC,
                                              float digitizationPeriod, float pedestal,
                                              TGraph &g, std::vector<std::pair<size_t,size_t> > &peaks)
{
  size_t nBins = waveform.size();
  size_t peakStartBin=0;
  size_t peakEndBin=0;
  for(size_t bin=0; bin<nBins; ++bin) 
  {
    g.SetPoint(bin,(startTDC+bin)*digitizationPeriod,waveform[bin]-pedestal);

    if(bin<2 || bin>=nBins-2) continue; //don't search for peaks here
    if(waveform[bin-1]<waveform[bin]) //rising edge
    {
      peakStartBin=bin;
      peakEndBin=bin;
    }
    if(waveform[bin-1]==waveform[bin]) //potentially a peak where consecutive ADC values are equal
    {
      peakEndBin=bin;
    }
    if(waveform[bin-1]>waveform[bin]) //falling edge
    {
      if(peakStartBin>0)  //found a peak
      {
        if(waveform[peakEndBin]-pedestal>_minADCdifference) peaks.emplace_back(peakStartBin,peakEndBin);
        peakStartBin=0;  //so that the loop has to wait for the next rising edge
      }
    }
  }
}

void MakeCrvRecoPulses::RangeFinderNarrow(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end)
{
  RangeFinder(waveform, peakStart, peakEnd, start, end);
  start=std::max(start,peakStart-4);
  end=std::min(end,peakEnd+4);
}

void MakeCrvRecoPulses::RangeFinder(const std::vector<unsigned int> &waveform, const size_t peakStart, const size_t peakEnd, size_t &start, size_t &end)
{
  for(size_t i=peakEnd; i<waveform.size()-1; ++i)
  {
    if(waveform[i+1]<=waveform[i]) end=i+1;
    else break;
  }
  for(size_t i=peakStart; i>0; --i)
  {
    if(waveform[i-1]<=waveform[i]) start=i-1;
    else break;
  }
}

bool MakeCrvRecoPulses::FailedFit(TFitResultPtr fr, int paramStart, int paramEnd)
{
  if(fr!=0) return true;
  if(!fr->IsValid()) return true;

  const double tolerance=0.01; //TODO: Need a user variable. Perhaps the limit condition can be extracted from the minimizer
  for(int i=paramStart; i<=paramEnd; ++i)
  {
    double v=fr->Parameter(i);
    double lower, upper;
    fr->ParameterBounds(i,lower,upper);
    if((v-lower)/(upper-lower)<tolerance) return true;
    if((upper-v)/(upper-lower)<tolerance) return true;
  }
  return false;
}

double MakeCrvRecoPulses::Chi2(TF1 &f, const TGraph &g)
{
  float chi2=0;
  double xmin,xmax;
  f.GetRange(xmin,xmax);
  for(int i=0; i<g.GetN(); ++i)
  {
    double x,y;
    g.GetPoint(i,x,y);
    if(x<xmin || x>xmax) continue;
    float fy=f.Eval(x);
    chi2+=(y-fy)*(y-fy)/fy;
  }
  return chi2;
}

void MakeCrvRecoPulses::NoFitOption(const std::vector<unsigned int> &waveform, float pedestal, 
                                    size_t peakStart, float &sum, size_t &pulseStart, size_t &pulseEnd)
{
  float maxADC=waveform[peakStart]-pedestal;
  float FWHMthreshold=maxADC/2.0;
  bool  foundPulseStart=false;
  bool  foundPulseEnd=false;
  sum=0;
  for(size_t i=peakStart+1; i<waveform.size(); ++i)
  {
    float ADC=waveform[i]-pedestal;
    sum+=ADC;
    if(ADC<FWHMthreshold && !foundPulseEnd) {pulseEnd=i; foundPulseEnd=true;}
    if(ADC<_minADCdifference) break;
  }
  if(!foundPulseEnd) pulseEnd=waveform.size()-1;
  for(size_t i=peakStart; ; --i)
  {
    float ADC=waveform[i]-pedestal;
    sum+=ADC;
    if(ADC<FWHMthreshold && !foundPulseStart) {pulseStart=i; foundPulseStart=true;}
    if(ADC<_minADCdifference || i==0) break;
  }
  if(!foundPulseStart) pulseStart=0;
}

void MakeCrvRecoPulses::SetWaveform(const std::vector<unsigned int> &waveform, 
                                    unsigned int startTDC, float digitizationPeriod, float pedestal, 
                                    float calibrationFactor, float calibrationFactorPulseHeight)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseBetas.clear();
  _pulseFitChi2s.clear();
  _PEs.clear();
  _PEsPulseHeight.clear();
  _LEtimes.clear();
  _failedFits.clear();
  _PEsNoFit.clear();
  _pulseTimesNoFit.clear();
  _pulseStart.clear();
  _pulseEnd.clear();

  //fill graph and find peaks
  std::vector<std::pair<size_t,size_t> > peaks;
  size_t nBins = waveform.size();
  TGraph g(nBins);
  FillGraphAndFindPeaks(waveform, startTDC, digitizationPeriod, pedestal, g, peaks);

  //loop through all peaks
  for(size_t ipeak=0; ipeak<peaks.size(); ++ipeak)
  {
    bool simpleFit=true;
    int  fittingNpeaks=1;

    size_t peakStartBin=peaks[ipeak].first;
    size_t peakEndBin=peaks[ipeak].second;
    double peakStartTime=(startTDC+peakStartBin)*digitizationPeriod;
    double peakEndTime=(startTDC+peakEndBin)*digitizationPeriod;
    double peakTime=0.5*(peakStartTime+peakEndTime);
    size_t peakStartBin2=0;
    size_t peakEndBin2=0;
    double peakStartTime2=0;
    double peakEndTime2=0;
    double peakTime2=0;

    size_t fitStartBin, fitEndBin;
    double fitStartTime, fitEndTime;

    //first try simple fit (=one Gumbel function)
    _f1.SetParameter(0, (waveform[peakStartBin]-pedestal)*TMath::E());
    _f1.SetParameter(1, peakTime);
    _f1.SetParameter(2, _defaultBeta);
    _f1.SetParLimits(0,(waveform[peakStartBin]-pedestal)*TMath::E()*_minPulseHeightRatio,(waveform[peakStartBin]-pedestal)*TMath::E()*_maxPulseHeightRatio);
    _f1.SetParLimits(1,peakStartTime-_maxTimeDifference,peakEndTime+_maxTimeDifference);
    _f1.SetParLimits(2, _minBeta, _maxBeta);

    RangeFinderNarrow(waveform, peakStartBin, peakEndBin, fitStartBin, fitEndBin);
    fitStartTime=(startTDC+fitStartBin)*digitizationPeriod;
    fitEndTime=(startTDC+fitEndBin)*digitizationPeriod;
    _f1.SetRange(fitStartTime,fitEndTime);

    //do the fit
    TFitResultPtr fr = g.Fit(&_f1,"NQSR");
//    float  pulseFitChi2 = _f1.GetChisquare()/_f1.GetNumberFitPoints();
    float  pulseFitChi2 = Chi2(_f1,g)/_f1.GetNumberFitPoints();
    double fitParam0 = fr->Parameter(0);
    double fitParam1 = fr->Parameter(1);
    double fitParam2 = fr->Parameter(2);
    double fitParam3 = 0;
    double fitParam4 = 0;
    double fitParam5 = 0;

    if((pulseFitChi2>_doubleGumbelThreshold || FailedFit(fr,0,2)) && _allowDoubleGumbel)
    {
      //try fit with two pulses (=two Gumbel function) (more time consuming)
      simpleFit=false;
      if(ipeak+1<peaks.size())
      {
        if(peaks[ipeak+1].first-peaks[ipeak].first<7) fittingNpeaks=2;  //two peaks close together
      }

      _f2.SetParameter(0, (waveform[peakStartBin]-pedestal)*TMath::E());
      _f2.SetParameter(1, peakTime);
      _f2.SetParameter(2, _defaultBeta);
      _f2.SetParameter(5, _defaultBeta);
      _f2.SetParLimits(0,(waveform[peakStartBin]-pedestal)*TMath::E()*_minPulseHeightRatio,(waveform[peakStartBin]-pedestal)*TMath::E()*_maxPulseHeightRatio);
      _f2.SetParLimits(1,peakStartTime-_maxTimeDifference,peakEndTime+_maxTimeDifference);
      _f2.SetParLimits(2, _minBeta, _maxBeta);
      _f2.SetParLimits(5, _minBeta, _maxBeta);
      if(fittingNpeaks==1)  //potentially merged double peaks or hidden second peak (most likely due to a reflected pulse)
      {
        _f2.SetParameter(3, 0.1*(waveform[peakStartBin]-pedestal)*TMath::E());
        _f2.SetParameter(4, peakTime+3*digitizationPeriod);  //TODO
        _f2.SetParLimits(3,0,(waveform[peakStartBin]-pedestal)*TMath::E()*_maxPulseHeightRatio);
        _f2.SetParLimits(4,peakEndTime+digitizationPeriod,peakEndTime+8*digitizationPeriod); //TODO

        size_t fitStartBin, fitEndBin;
        RangeFinder(waveform, peakStartBin, peakEndBin, fitStartBin, fitEndBin);
        double fitStartTime=(startTDC+fitStartBin)*digitizationPeriod;
        double fitEndTime=(startTDC+fitEndBin)*digitizationPeriod;
        _f2.SetRange(fitStartTime,fitEndTime);
      }
      else //separated double peaks (most likely due to a reflected pulse)
      {
        peakStartBin2=peaks[ipeak+1].first;
        peakEndBin2=peaks[ipeak+1].second;
        peakStartTime2=(startTDC+peakStartBin2)*digitizationPeriod;
        peakEndTime2=(startTDC+peakEndBin2)*digitizationPeriod;
        peakTime2=0.5*(peakStartTime2+peakEndTime2);

        _f2.SetParameter(3, (waveform[peakStartBin2]-pedestal)*TMath::E());
        _f2.SetParameter(4, peakTime2);
        _f2.SetParLimits(3,(waveform[peakStartBin2]-pedestal)*TMath::E()*_minPulseHeightRatio,(waveform[peakStartBin2]-pedestal)*TMath::E()*_maxPulseHeightRatio);
        _f2.SetParLimits(4,peakStartTime2-_maxTimeDifference,peakEndTime2+_maxTimeDifference);

        size_t fitStartBin, fitEndBin;
        RangeFinder(waveform, peakStartBin, peakEndBin, fitStartBin, fitEndBin);
        double fitStartTime=(startTDC+fitStartBin)*digitizationPeriod;
        RangeFinder(waveform, peakStartBin2, peakEndBin2, fitStartBin, fitEndBin);
        double fitEndTime=(startTDC+fitEndBin)*digitizationPeriod;
        _f2.SetRange(fitStartTime,fitEndTime);

        ++ipeak;  //this loop already looked at the next peaks, so we need to skip it in the next loop
      }

      //do the fit
      fr = g.Fit(&_f2,"NQSR");
//      pulseFitChi2 = _f2.GetChisquare()/_f2.GetNumberFitPoints();
      pulseFitChi2 = Chi2(_f2,g)/_f2.GetNumberFitPoints();
      fitParam0 = fr->Parameter(0);
      fitParam1 = fr->Parameter(1);
      fitParam2 = fr->Parameter(2);
      fitParam3 = fr->Parameter(3);
      fitParam4 = fr->Parameter(4);
      fitParam5 = fr->Parameter(5);
    } //fit with two Gumbel functions

    //collect fit information for the first peak
    float  PEs          = fitParam0*fitParam2 / calibrationFactor;
    double pulseTime    = fitParam1;
    float  pulseHeight  = fitParam0/TMath::E();
    float  pulseBeta    = fitParam2;
    double LEtime       = 0;
    bool   failedFit    = FailedFit(fr,0,2);

    if(failedFit)
    {
      PEs          = (waveform[peakStartBin]-pedestal)*TMath::E() * _defaultBeta / calibrationFactor;
      pulseTime    = peakTime;
      pulseHeight  = waveform[peakStartBin]-pedestal;
      pulseBeta    = _defaultBeta;
      pulseFitChi2 = NAN;
      LEtime       = peakTime-_LEtimeFactor*_defaultBeta;  //50% pulse height is reached at -0.985*beta before the peak
      failedFit    = true;
    }
    else LEtime    = pulseTime-_LEtimeFactor*pulseBeta;  //50% pulse height is reached at -0.985*beta before the peak

    float PEsPulseHeight = pulseHeight / calibrationFactorPulseHeight;

    _pulseTimes.push_back(pulseTime);
    _pulseHeights.push_back(pulseHeight);
    _pulseBetas.push_back(pulseBeta);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _PEs.push_back(PEs);
    _PEsPulseHeight.push_back(PEsPulseHeight);
    _LEtimes.push_back(LEtime);
    _failedFits.push_back(failedFit);

    float  sum;
    size_t pulseStartBin;
    size_t pulseEndBin;
    NoFitOption(waveform, pedestal, peakStartBin, sum, pulseStartBin, pulseEndBin);
    _PEsNoFit.push_back(sum*digitizationPeriod/calibrationFactor);
    _pulseTimesNoFit.push_back(peakTime);
    _pulseStart.push_back((startTDC+pulseStartBin)*digitizationPeriod);
    _pulseEnd.push_back((startTDC+pulseEndBin)*digitizationPeriod);

//std::cout<<_PEs.back()<<"  "<<_PEsNoFit.back()<<"          "<<_pulseStart.back()<<"  "<<_pulseTimes.back()<<" / "<<_pulseTimesNoFit.back()<<"  "<<_pulseEnd.back()<<std::endl;

    if(simpleFit) continue;  //no second peak

    //collect fit information for the second peak
    float  PEs_2          = fitParam3*fitParam5 / calibrationFactor;
    double pulseTime_2    = fitParam4;
    float  pulseHeight_2  = fitParam3/TMath::E();
    float  pulseBeta_2    = fitParam5;
    double LEtime_2       = pulseTime_2-_LEtimeFactor*pulseBeta_2;  //50% pulse height is reached at -0.985*beta before the peak
    bool   failedFit_2    = FailedFit(fr,3,5);
    float  PEsPulseHeight_2 = pulseHeight_2 / calibrationFactorPulseHeight;

    if(fittingNpeaks==1)
    {
      if(failedFit) continue; //the main peak didn't work, no need to continue with a hidden second peak
      if(failedFit_2) continue; //there was no hidden second peak
      if(!failedFit_2 && PEs_2<1.0) continue; //there was no hidden second peak
    }
    if(fittingNpeaks==2)
    {
      if(failedFit_2) {continue; --ipeak;} //repeat this peak
    }

    _pulseTimes.push_back(pulseTime_2);
    _pulseHeights.push_back(pulseHeight_2);
    _pulseBetas.push_back(pulseBeta_2);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _PEs.push_back(PEs_2);
    _PEsPulseHeight.push_back(PEsPulseHeight_2);
    _LEtimes.push_back(LEtime_2);
    _failedFits.push_back(failedFit_2);

/*
    _PEs.back()+=PEs_2;
    _PEsPulseHeight.back()+=PEsPulseHeight_2;
*/

    if(fittingNpeaks==2)
    {
      NoFitOption(waveform, pedestal, peakStartBin2, sum, pulseStartBin, pulseEndBin);
      _PEsNoFit.push_back(sum*digitizationPeriod/calibrationFactor);
      _pulseTimesNoFit.push_back(peakTime2);
      _pulseStart.push_back((startTDC+pulseStartBin)*digitizationPeriod);
      _pulseEnd.push_back((startTDC+pulseEndBin)*digitizationPeriod);
    }
    else
    {
      _PEsNoFit.push_back(0);
      _pulseTimesNoFit.push_back(0);
      _pulseStart.push_back(0);
      _pulseEnd.push_back(0);
    }
    
  }

}

}

