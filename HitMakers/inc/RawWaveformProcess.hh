///////////////////////////////////////////////////////////////////////////////
// 2016-01-21 G. Pezzullo
// -------------------
// resolve waveform and creates RecoCaloDigi
///////////////////////////////////////////////////////////////////////////////

#ifndef RawWaveformProcess_HH
#define RawWaveformProcess_HH

#include "HitMakers/inc/WaveformProcess.hh"


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
  struct Hist_raw2 {
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
    TH1F*  _hFitM;
      
      
    TH1F*  _debugWf[20];
    TH2F*  _hEdepAmp;
  };

  class RawWaveformProcess : public WaveformProcess {
     
    

  protected:
    //external parameters to set at the initialization stage
    int         _debugLevel;
    int         _acquisitionEndTime;
    int         _digiSampling;

    //internal variables used in the methods
    double      _ADCToMeV;
    double      _ADCPeak2MeV;
    double      _psd;       //pulse shape discriminator 
    int         _nMax;      //number of maxima found in the pulse
    double      _max [20];  //array of the maxima
    double      _eDep[20];  //pulse from the peak
    double      _time[20];  //time associated to the peak

    double      _amplitude, _charge;

    Hist_raw2   _hist;
    TTree*      _tree;

    //    int         _debugHistIndex;

    Int_t       _counter, _nParticles, _hitCounter, _pdgId[1000], _nPeaks;
    Float_t     _timeWf, _mcTime, _mcMeanTime, _mcEDep, _psdWf,  _Chi2Time, _nDof, _fitM, _fitQ, _amp;

//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  public:
#ifndef __GCCXML__
    explicit    RawWaveformProcess(fhicl::ParameterSet const&);
#endif/*__GCCXML__*/
    virtual    ~RawWaveformProcess();
    
    virtual void      processWaveform(double ADCToMeV, CaloDigi CaloHit, RecoCaloDigiCollection &RecoCaloHits);


    virtual void      fillDiag(CaloDigiMC*DigiMC, RecoCaloDigiCollection *RecoCaloHits);

    virtual void      book();

  };
}
#endif
