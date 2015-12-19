///////////////////////////////////////////////////////////////////////////////
// 2015-04-08 G. Pezzullo
// -------------------
// resolve waveform and creates RecoCaloDigi
///////////////////////////////////////////////////////////////////////////////

#ifndef FitWaveformProcess_HH
#define FitWaveformProcess_HH

#include "HitMakers/inc/WaveformProcess.hh"

// #include "TH1F.h"
// #include "TF1.h"
// #include "TFile.h"
// #include "TROOT.h"
// #include "TFolder.h"
// #include "TDirectory.h"


#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/
#include <cstddef>
#include <vector>
#include "Rtypes.h"

namespace art { class TFileDirectory; }

class TH1F;
class TH2F;
class TTree;
class TFile;
class TDirectory;
class TF1;
class TFolder;

namespace mu2e {

  struct Hist_t {
    TH1F*  _hEdep;
    TH1F*  _hAmplitude;
    TH1F*  _hTime;

    TH2F*  _hTimeVsE;
    TH2F*  _hTimeVsAmp;
    TH2F*  _hAmpVsE;
    TH2F*  _hPSDVsAmp;
    
    TH1F*  _hDt;
    TH1F*  _hPSD;
    TH1F*  _hChi2Time;
    TH1F*  _hNDofTime;
    TH1F*  _hDtBinTime;
    TH1F*  _hFitEta;
    TH1F*  _hFitPeak;
    TH1F*  _hFitSigma;
    TH1F*  _hFitNorm;
    TH1F*  _debugWf[20];
    TH2F*  _hEdepAmp;
    
    

  };

  class FitWaveformProcess : public WaveformProcess {
     
    

  protected:
    //external parameters to set at the initialization stage
    int         _debugLevel;
    double      _acquisitionLenght;
    int         _digiSampling;
    int         _wave_point_error;
    double      _psdThreshold;
    double      _timeFraction;
    
    //internal variables used in the methods
    double      _ADCToMeV;
    double      _psd;       //pulse shape discriminator 
    int         _nMax;      //number of maxima found in the pulse
    double      _max [20];   //array of the maxima
    double      _eDep[2];
    double      _time[2];
    double      _timeChi2[2];
    TH1F*       _wave;      //input waveform
    TH1F*       _refPulse;  //pulse used as reference
    TF1*        _flogn;
    TF1*        _fdlogn;

    Hist_t      _hist;
    TTree*      _tree;

    int         _debugHistIndex;

    Int_t       _counter, _nParticles, _hitCounter, _pdgId[1000];
    Float_t     _timeWf, _mcTime, _mcEDep, _psdWf,  _Chi2Time, _nDof, _fitEta, _fitNorm, _fitSigma, _fitPeak,_fitDlogn, _amp, _charge;

//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  public:
#ifndef __GCCXML__
    explicit    FitWaveformProcess(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
    virtual    ~FitWaveformProcess();
    
    //check for presence of pile up, fills _psd
    void        psdFromChi2();
    
    void        psdFromAmplitude();
    
    //search maxima and minima, filling the protected variables: _nMax, _nMin, _max[], _min[]
    void        findMaxima();
    
    //perform a fit for separating the contributes in the pulse
    void        separeteWaveform();

    //integrate the pulse
    void        calculateEnergy (double &Edep);
    
    //perform a fit only on the leading edge and use digital-constant-fraction for time extraction
    void        calculateTime   (double WaveMax, double &Time, double &Chi2);
    

    //asymmetric log-Normal from Grupen
    static double     logn(double *x, double *par);

    //sum of two log-normal for identifying two pulses
    static double     dlogn(double *x, double *par);



    virtual void      processWaveform(double ADCToMeV, CaloDigi CaloHit, RecoCaloDigiCollection &RecoCaloHits);


    virtual void      fillDiag(CaloDigiMC*DigiMC, RecoCaloDigiCollection *RecoCaloHits);

    virtual void      book();

  };
}
#endif
