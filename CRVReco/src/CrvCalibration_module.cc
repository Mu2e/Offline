//
// A module that finds the calibration constants for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/DataProducts/inc/CRVId.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TSpectrum.h>

namespace mu2e
{
  class CrvCalibration : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{Name("crvRecoPulsesModuleLabel"), Comment("module label of the input CrvRecoPulses")};
      fhicl::Atom<bool>        useNZS{Name("useNZS"), Comment("use NZS data"), false};
      fhicl::Atom<int>         histBinsPulseArea{Name("histBinsPulseArea"), Comment("pulseArea histogram bins"), 300};
      fhicl::Atom<int>         histBinsPulseHeight{Name("histBinsPulseHeight"), Comment("pulseHeight histogram bins"), 300};
      fhicl::Atom<double>      histMaxPulseArea{Name("histMaxPulseArea"), Comment("end range of pulseArea histogram"), 3000.0};
      fhicl::Atom<double>      histMaxPulseHeight{Name("histMaxPulseHeight"), Comment("end range of pulseArea histogram"), 150.0};
      fhicl::Atom<double>      fitRangeStart{Name("fitRangeStart"), Comment("low end of the 1PE fit range as fraction of peak"), 0.8};
      fhicl::Atom<double>      fitRangeEnd{Name("fitRangeEnd"), Comment("high end of the 1PE fit range as fraction of peak"), 1.2};
      fhicl::Atom<bool>        use2PE1PEdifference{Name("use2PE1PEdifference"), Comment("use the 2PE-1PE difference as calib constant"), false}; //removes impact of pedestal offsets, but returns unstable results
      fhicl::Atom<double>      fitRangeStart2PE{Name("fitRangeStart2PE"), Comment("low end of the 2PE fit range as fraction of peak"), 0.9};
      fhicl::Atom<double>      fitRangeEnd2PE{Name("fitRangeEnd2PE"), Comment("high end of the 2PE fit range as fraction of peak"), 1.1};
      fhicl::Atom<double>      minPeakPulseArea{Name("minPeakPulseArea"), Comment("minimum accepted 1PE peak for pulseArea histogram"), 280.0};
      fhicl::Atom<double>      minPeakPulseHeight{Name("minPeakPulseHeight"), Comment("minimum accepted 1PE peak for pulseHeight histogram"), 10.0};
      fhicl::Atom<int>         minHistEntries{Name("minHistEntries"), Comment("minimum number of entries required for a fit"), 100};
      fhicl::Atom<int>         spectrumNPeaks{Name("spectrumNPeaks"), Comment("maximum number of peaks searched by TSpectrum"), 100};
      fhicl::Atom<double>      spectrumPeakSigma{Name("spectrumPeakSigma"), Comment("TSpectrum search parameter sigma"), 4.0};
      fhicl::Atom<double>      spectrumPeakThreshold{Name("spectrumPeakThreshold"), Comment("TSpectrum search parameter threshold"), 0.01};
      fhicl::Atom<double>      maxFitDifferencePulseArea{Name("maxFitDifferencePulseArea"),
                                   Comment("maximum difference between the TSpectrum peak and the fitted peak (for pulse areas). Indicates fit problems."), 80.0};
      fhicl::Atom<double>      maxFitDifferencePulseHeight{Name("maxFitDifferencePulseHeight"),
                                   Comment("maximum difference between the TSpectrum peak and the fitted peak (for pulse heights). Indicates fit problems."), 3.0};
      fhicl::Atom<std::string> tmpDBfileName{Name("tmpDBfileName"), Comment("name of the tmp. DB file name for the pedestals")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvCalibration(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);
    void endJob();
    bool Find1PEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &OnePEpeak, double minPeak, double maxFitDifference);
    bool Find2PEpeak(TH1F *hist, TF1 &function2PE, double OnePEpeak, double &TwoPEpeak);

    private:
    std::string        _crvRecoPulsesModuleLabel;
    bool               _useNZS;
    int                _histBinsPulseArea, _histBinsPulseHeight;
    double             _histMaxPulseArea, _histMaxPulseHeight;
    double             _fitRangeStart, _fitRangeEnd;
    bool               _use2PE1PEdifference;
    double             _fitRangeStart2PE, _fitRangeEnd2PE;
    int                _minHistEntries;
    double             _minPeakPulseArea, _minPeakPulseHeight;
    int                _spectrumNPeaks;
    double             _spectrumPeakSigma;
    double             _spectrumPeakThreshold;
    double             _maxFitDifferencePulseArea, _maxFitDifferencePulseHeight;
    std::string        _tmpDBfileName;
    std::vector<TH1F*> _calibHistsPulseArea;
    std::vector<TH1F*> _calibHistsPulseHeight;

    ProditionsHandle<CRVCalib> _calib_h;

    std::vector<double> _pedestals;
    bool                _first{true};
  };


  CrvCalibration::CrvCalibration(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _useNZS(conf().useNZS()),
    _histBinsPulseArea(conf().histBinsPulseArea()),
    _histBinsPulseHeight(conf().histBinsPulseHeight()),
    _histMaxPulseArea(conf().histMaxPulseArea()),
    _histMaxPulseHeight(conf().histMaxPulseHeight()),
    _fitRangeStart(conf().fitRangeStart()),
    _fitRangeEnd(conf().fitRangeEnd()),
    _use2PE1PEdifference(conf().use2PE1PEdifference()),
    _fitRangeStart2PE(conf().fitRangeStart2PE()),
    _fitRangeEnd2PE(conf().fitRangeEnd2PE()),
    _minHistEntries(conf().minHistEntries()),
    _minPeakPulseArea(conf().minPeakPulseArea()),
    _minPeakPulseHeight(conf().minPeakPulseHeight()),
    _spectrumNPeaks(conf().spectrumNPeaks()),
    _spectrumPeakSigma(conf().spectrumPeakSigma()),
    _spectrumPeakThreshold(conf().spectrumPeakThreshold()),
    _maxFitDifferencePulseArea(conf().maxFitDifferencePulseArea()),
    _maxFitDifferencePulseHeight(conf().maxFitDifferencePulseHeight()),
    _tmpDBfileName(conf().tmpDBfileName())
  {
  }

  void CrvCalibration::beginRun(const art::Run&)
  {
    if(_calibHistsPulseArea.size()>0) return;  //don't initialize again for additional runs

    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    _calibHistsPulseArea.reserve(counters.size()*CRVId::nChanPerBar);
    _calibHistsPulseHeight.reserve(counters.size()*CRVId::nChanPerBar);
    _pedestals.resize(counters.size()*CRVId::nChanPerBar);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
    {
      for(size_t SiPM=0; SiPM<CRVId::nChanPerBar; ++SiPM)
      {
        //produce histograms also for non-existing channels to get a continuously running index
        size_t channelIndex=barIndex*CRVId::nChanPerBar + SiPM;
        _calibHistsPulseArea.emplace_back(tfs->make<TH1F>(Form("crvCalibrationHistPulseArea_%lu",channelIndex),
                                                 Form("crvCalibrationHistPulseArea_%lu",channelIndex),
                                                 _histBinsPulseArea,0,_histMaxPulseArea));
        _calibHistsPulseHeight.emplace_back(tfs->make<TH1F>(Form("crvCalibrationHistPulseHeight_%lu",channelIndex),
                                                 Form("crvCalibrationHistPulseHeight_%lu",channelIndex),
                                                 _histBinsPulseHeight,0,_histMaxPulseHeight));
        _pedestals[channelIndex]=0;
      }
    }
  }

  void CrvCalibration::endJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    TTree *treePedestal = tfs->make<TTree>("crvPedestals","crvPedestals");
    size_t channel;
    double pedestal;
    treePedestal->Branch("channel", &channel);
    treePedestal->Branch("pedestal", &pedestal);

    TF1 function("calibPeak1PE","gaus");
    TF1 function2PE("calibPeak2PE","gaus");
    TSpectrum spectrum(_spectrumNPeaks);

    std::ofstream outputFile;
    outputFile.open(_tmpDBfileName);
    if(!outputFile.is_open()) throw cet::exception("CRVCALIB") << "Couldn't create output txt file " << _tmpDBfileName << ".";
    outputFile<<"TABLE CRVSiPM"<<std::endl;
    outputFile<<"#channel,pedestal,calibPulseHeight,calibPulseArea"<<std::endl;

    for(channel=0; channel<_pedestals.size(); ++channel)
    {
      TH1F *hist;
      double calibValue[2]={-1,-1};
      for(int i=0; i<2; ++i) //loop over hisograms with pulse areas and pulse heights
      {
        if(i==1) hist=_calibHistsPulseArea.at(channel);
        else hist=_calibHistsPulseHeight.at(channel);

        double OnePEpeak=-1;
        if(!Find1PEpeak(hist, spectrum, function, OnePEpeak, (i==0?_minPeakPulseHeight:_minPeakPulseArea), (i==0?_maxFitDifferencePulseHeight:_maxFitDifferencePulseArea))) continue;
        if(_use2PE1PEdifference)
        {
          double TwoPEpeak=-1;
          if(!Find2PEpeak(hist, function2PE, OnePEpeak, TwoPEpeak)) continue;
          calibValue[i]=TwoPEpeak-OnePEpeak;
        }
        else calibValue[i]=OnePEpeak;
      }

      pedestal=_pedestals.at(channel);
      outputFile<<channel<<","<<pedestal<<","<<calibValue[0]<<","<<calibValue[1]<<std::endl;  //write to DB text file
      treePedestal->Fill(); //fill tree
    }

    outputFile.close();
  }

  void CrvCalibration::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,(_useNZS?"NZS":""),crvRecoPulseCollection)) return;

    //find pedestals and time offsets from first event
    //need to assume that this is only used for calibration runs where both values stay constant over the entire run
    if(_first)
    {
      _first=false;
      auto const& calib = _calib_h.get(event.id());
      for(size_t channelIndex=0; channelIndex<_pedestals.size(); ++channelIndex)
      {
        _pedestals[channelIndex] = calib.pedestal(channelIndex);
      }
    }

    for(auto iter=crvRecoPulseCollection->begin(); iter!=crvRecoPulseCollection->end(); ++iter)
    {
      //check for any error flags, but ignore calibration flags
      auto flags = iter->GetRecoPulseFlags();
      flags.set(CrvRecoPulseFlagEnums::noCalibConstPulseArea,0);
      flags.set(CrvRecoPulseFlagEnums::noCalibConstPulseHeight,0);
      if(!flags.none()) continue;

      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      int channelIndex=barIndex*CRVId::nChanPerBar+SiPM;

      _calibHistsPulseArea.at(channelIndex)->Fill(iter->GetPulseBeta()*iter->GetPulseHeight()*TMath::E());
      _calibHistsPulseHeight.at(channelIndex)->Fill(iter->GetPulseHeight());
    }
  }

  bool CrvCalibration::Find1PEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &OnePEpeak, double minPeak, double maxFitDifference)
  {
    if(hist->GetEntries()<_minHistEntries) return false; //not enough data

    int nPeaks = spectrum.Search(hist,_spectrumPeakSigma,"nodraw",_spectrumPeakThreshold);
    if(nPeaks<=0) return false;

    //peaks are returned sorted by Y
    //from our long-time experience:
    //-1PE peak is either the highest peak or second highest peak
    //-if the 1PE peak is the second highest, then the highest peak comes from the baseline
    //-the peak from the baseline is always below the minPeak threshold, while the 1PE peak is not
    //-the minPeak threshold may have to be adjusted for non-standard bias voltages
    double *peaksX = spectrum.GetPositionX();
    double x=peaksX[0];
    if(x<minPeak)
    {
      if(nPeaks==1) return false;
      x=peaksX[1];
      if(x<minPeak) return false;
    }

    if(hist->FindBin(x*_fitRangeStart)==hist->FindBin(x*_fitRangeEnd)) return false; //fit range start/end are in the same bin
    function.SetRange(x*_fitRangeStart,x*_fitRangeEnd);
    function.SetParameter(1,x);
    hist->Fit(&function, "QR");
    OnePEpeak = function.GetParameter(1);

    if(fabs(OnePEpeak-x)>maxFitDifference) return false;
    if(OnePEpeak<minPeak) return false;
    return true;
  }

  bool CrvCalibration::Find2PEpeak(TH1F *hist, TF1 &function2PE, double OnePEpeak, double &TwoPEpeak)
  {
    if(hist->GetEntries()<_minHistEntries) return false; //not enough data

    if(OnePEpeak<=0) return false;  //no valid 1PE peak

    hist->GetXaxis()->SetRangeUser(2.0*OnePEpeak*_fitRangeStart2PE,2.0*OnePEpeak*_fitRangeEnd2PE);
    float x = hist->GetBinCenter(hist->GetMaximumBin());
    hist->GetXaxis()->SetRange(0,0); //UnZoom()

    function2PE.SetRange(x*_fitRangeStart2PE,x*_fitRangeEnd2PE);
    function2PE.SetParameter(1,x);
    hist->Fit(&function2PE, "QR");
    TwoPEpeak = function2PE.GetParameter(1);

    if(TwoPEpeak<x*_fitRangeStart2PE || TwoPEpeak>x*_fitRangeEnd2PE) return false;
    return true;
  }

} // end namespace mu2e

using mu2e::CrvCalibration;
DEFINE_ART_MODULE(CrvCalibration)
