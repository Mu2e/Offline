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
      fhicl::Atom<int>         histBinsPulseArea{Name("histBinsPulseArea"), Comment("pulseArea histogram bins"), 300};
      fhicl::Atom<int>         histBinsPulseHeight{Name("histBinsPulseHeight"), Comment("pulseHeight histogram bins"), 300};
      fhicl::Atom<double>      histMaxPulseArea{Name("histMaxPulseArea"), Comment("end range of pulseArea histogram"), 3000.0};
      fhicl::Atom<double>      histMaxPulseHeight{Name("histMaxPulseHeight"), Comment("end range of pulseArea histogram"), 150.0};
      fhicl::Atom<double>      fitRangeStart{Name("fitRangeStart"), Comment("low end of the 1PE fit range as fraction of peak"), 0.8};
      fhicl::Atom<double>      fitRangeEnd{Name("fitRangeEnd"), Comment("high end of the 1PE fit range as fraction of peak"), 1.2};
      fhicl::Atom<double>      minPeakPulseArea{Name("minPeakPulseArea"), Comment("minimum accepted SPE peak for pulseArea histogram"), 250.0};
      fhicl::Atom<double>      minPeakPulseHeight{Name("minPeakPulseHeight"), Comment("minimum accepted SPE peak for pulseHeight histogram"), 10.0};
      fhicl::Atom<int>         minHistEntries{Name("minHistEntries"), Comment("minimum number of entries required for a fit"), 100};
      fhicl::Atom<int>         spectrumNPeaks{Name("spectrumNPeaks"), Comment("maximum number of peaks searched by TSpectrum"), 100};
      fhicl::Atom<double>      spectrumPeakSigma{Name("spectrumPeakSigma"), Comment("TSpectrum search parameter sigma"), 4.0};
      fhicl::Atom<double>      spectrumPeakThreshold{Name("spectrumPeakThreshold"), Comment("TSpectrum search parameter threshold"), 0.01};
      //fhicl::Atom<double>      peakRatioTolerance{Name("peakRatioTolerance"), Comment("allowed deviation of the ratio between 1PE peak and 2PE peak from 2.0"), 0.2};
      fhicl::Atom<std::string> tmpDBfileName{Name("tmpDBfileName"), Comment("name of the tmp. DB file name for the pedestals")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvCalibration(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);
    void endJob();
//    template<size_t N>
//    bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, std::array<TF1,N> &functions, double &SPEpeak, double minPeak);
    bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &SPEpeak, double minPeak);

    private:
    std::string        _crvRecoPulsesModuleLabel;
    int                _histBinsPulseArea, _histBinsPulseHeight;
    double             _histMaxPulseArea, _histMaxPulseHeight;
    double             _fitRangeStart, _fitRangeEnd;
    int                _minHistEntries;
    double             _minPeakPulseArea, _minPeakPulseHeight;
    int                _spectrumNPeaks;
    double             _spectrumPeakSigma;
    double             _spectrumPeakThreshold;
//    double             _peakRatioTolerance;
    std::string        _tmpDBfileName;
    std::vector<TH1F*> _calibHistsPulseArea;
    std::vector<TH1F*> _calibHistsPulseHeight;

    ProditionsHandle<CRVCalib> _calib_h;

    std::vector<double> _pedestals;
    std::vector<double> _timeOffsets;

    std::pair<int,int>  _firstRunSubrun;
    std::pair<int,int>  _lastRunSubrun;
  };


  CrvCalibration::CrvCalibration(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _histBinsPulseArea(conf().histBinsPulseArea()),
    _histBinsPulseHeight(conf().histBinsPulseHeight()),
    _histMaxPulseArea(conf().histMaxPulseArea()),
    _histMaxPulseHeight(conf().histMaxPulseHeight()),
    _fitRangeStart(conf().fitRangeStart()),
    _fitRangeEnd(conf().fitRangeEnd()),
    _minHistEntries(conf().minHistEntries()),
    _minPeakPulseArea(conf().minPeakPulseArea()),
    _minPeakPulseHeight(conf().minPeakPulseHeight()),
    _spectrumNPeaks(conf().spectrumNPeaks()),
    _spectrumPeakSigma(conf().spectrumPeakSigma()),
    _spectrumPeakThreshold(conf().spectrumPeakThreshold()),
//    _peakRatioTolerance(conf().peakRatioTolerance()),
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
    _timeOffsets.resize(counters.size()*CRVId::nChanPerBar);

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
        _timeOffsets[channelIndex]=0;
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

//    std::array functions={TF1("calibPeak1","gaus"), TF1("calibPeak2","gaus"), TF1("calibPeak3","gaus")}; //only need to fit three peaks
    TF1 function("calibPeak","gaus");
    TSpectrum spectrum(_spectrumNPeaks); //any value of 3 or less results in a "Peak buffer full" warning.

    std::ofstream outputFile;
    outputFile.open(_tmpDBfileName);
    outputFile<<"TABLE CRVSiPM "<<_firstRunSubrun.first<<":"<<_firstRunSubrun.second<<"-"<<_lastRunSubrun.first<<":"<<_lastRunSubrun.second<<std::endl;
    outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

    for(channel=0; channel<_pedestals.size(); ++channel)
    {
      TH1F *hist;
      double calibValue[2];
      for(int i=0; i<2; ++i) //loop over hisograms with pulse areas and pulse heights
      {
        if(i==1) hist=_calibHistsPulseArea.at(channel);
        else hist=_calibHistsPulseHeight.at(channel);

        double SPEpeak=-1;
//        if(!FindSPEpeak(hist, spectrum, functions, SPEpeak, (i==0?_minPeakPulseHeight:_minPeakPulseArea)))
        if(!FindSPEpeak(hist, spectrum, function, SPEpeak, (i==0?_minPeakPulseHeight:_minPeakPulseArea)))
        {
          calibValue[i]=-1;
          continue;
        }
        calibValue[i]=SPEpeak;
      }

      pedestal=_pedestals.at(channel);
      outputFile<<channel<<","<<pedestal<<","<<calibValue[0]<<","<<calibValue[1]<<std::endl;  //write to DB text file
      treePedestal->Fill(); //fill tree
    }

    outputFile<<std::endl;

    //time offsets
    TTree *treeTimeOffset = tfs->make<TTree>("crvTimeOffsets","crvTimeOffsets");
    double offset;
    treeTimeOffset->Branch("channel", &channel);
    treeTimeOffset->Branch("timeOffset", &offset);

    outputFile<<"TABLE CRVTime"<<std::endl;
    outputFile<<"#channel, timeOffset"<<std::endl;
    for(channel=0; channel<_timeOffsets.size(); ++channel)
    {
      offset=_timeOffsets.at(channel);
      outputFile<<channel<<","<<offset<<std::endl;  //write to DB text file
      treeTimeOffset->Fill(); //fill tree
    }

    outputFile.close();
  }

  void CrvCalibration::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"NZS",crvRecoPulseCollection)) return;

    //find pedestals and time offsets from first event
    //need to assume that this is only used for calibration runs where both values stay constant over the entire run
    static bool first=true;
    if(first)
    {
      first=false;
      auto const& calib = _calib_h.get(event.id());
      for(size_t channelIndex=0; channelIndex<_pedestals.size(); ++channelIndex)
      {
        _pedestals[channelIndex] = calib.pedestal(channelIndex);
        _timeOffsets[channelIndex] = calib.timeOffset(channelIndex);
      }

      _firstRunSubrun=std::pair<int,int>(event.run(),event.subRun());
    }
    _lastRunSubrun=std::pair<int,int>(event.run(),event.subRun());

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

  bool CrvCalibration::FindSPEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &SPEpeak, double minPeak)
  {
    if(hist->GetEntries()<_minHistEntries) return false; //not enough data

    int nPeaks = spectrum.Search(hist,_spectrumPeakSigma,"nodraw",_spectrumPeakThreshold);
    if(nPeaks<=0) return false;

    //peaks are returned sorted by Y
    //from our long-time experience:
    //-SPE peak is either the highest peak or second highest peak
    //-if the SPE peak is the second highest, then the highest peak comes from the baseline
    //-the peak from the baseline is always below the minPeak threshold, while the SPE peak is not
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
    SPEpeak = function.GetParameter(1);

    return true;
  }

/*
  template<size_t N>
  bool CrvCalibration::FindSPEpeak(TH1F *hist, TSpectrum &spectrum, std::array<TF1,N> &functions, double &SPEpeak, double minPeak)
  {
    if(hist->GetEntries()<_minHistEntries) return false; //not enough data

    size_t nPeaks = spectrum.Search(hist,_spectrumPeakSigma,"nodraw",_spectrumPeakThreshold);
    if(nPeaks==0) return false;

    //peaks are returned sorted by Y
    double *peaksX = spectrum.GetPositionX();
    std::vector<double> fittedPeaks;
    for(size_t iPeak=0; iPeak<nPeaks && iPeak<functions.size(); ++iPeak)
    {
      double x=peaksX[iPeak];
      if(hist->FindBin(x*_fitRangeStart)==hist->FindBin(x*_fitRangeEnd)) continue; //fit range start/end are in the same bin
      functions[iPeak].SetRange(x*_fitRangeStart,x*_fitRangeEnd);
      functions[iPeak].SetParameter(1,x);
      hist->Fit(&functions[iPeak], "QR+");
      fittedPeaks.emplace_back(functions[iPeak].GetParameter(1));
    }
    if(fittedPeaks.size()==0) return false;

    int peakToUse=-1;
    //only need to test two highest peaks (=first two entries in vector)
    //one of the peaks could be due to baseline fluctuations
    for(size_t iPeak=0; iPeak<fittedPeaks.size() && iPeak<2; ++iPeak)
    for(size_t jPeak=iPeak+1; jPeak<fittedPeaks.size(); ++jPeak)
    {
      if(fabs(fittedPeaks.at(jPeak)/fittedPeaks.at(iPeak)-2.0)<_peakRatioTolerance) peakToUse=iPeak;
    }
    if(peakToUse==-1) return false;
    if(fittedPeaks.at(peakToUse)<minPeak) return false;

    SPEpeak = fittedPeaks.at(peakToUse);

    return true;
  }
*/

} // end namespace mu2e

using mu2e::CrvCalibration;
DEFINE_ART_MODULE(CrvCalibration)
