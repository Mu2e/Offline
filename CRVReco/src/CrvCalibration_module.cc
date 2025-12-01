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
      fhicl::Atom<int>         histBinsPulseArea{Name("histBinsPulseArea"), Comment("pulseArea histogram bins"), 150};
      fhicl::Atom<int>         histBinsPulseHeight{Name("histBinsPulseHeight"), Comment("pulseHeight histogram bins"), 150};
      fhicl::Atom<double>      histMaxPulseArea{Name("histMaxPulseArea"), Comment("end range of pulseArea histogram"), 3000.0};
      fhicl::Atom<double>      histMaxPulseHeight{Name("histMaxPulseHeight"), Comment("end range of pulseArea histogram"), 150.0};
      fhicl::Atom<std::string> tmpDBfileName{Name("tmpDBfileName"), Comment("name of the tmp. DB file name for the pedestals")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvCalibration(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);
    void endJob();

    private:
    std::string        _crvRecoPulsesModuleLabel;
    int                _histBinsPulseArea, _histBinsPulseHeight;
    double             _histMaxPulseArea, _histMaxPulseHeight;
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

    TF1 funcCalib("f0", "gaus");

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

        if(hist->GetEntries()<100) //not enough data
        {
          calibValue[i]=-1;
          continue;
        }

/*
        int n=hist->GetNbinsX();
        double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
        if(overflow/((double)hist->GetEntries())>0.1) //too much underflow/overflow. something may be wrong.
        {
          calibValue[i]=-1;
          continue;
        }
*/

        int maxbinCalib = hist->GetMaximumBin();
        double peakCalib = hist->GetBinCenter(maxbinCalib);
//FIXME        funcCalib.SetRange(peakCalib*0.8,peakCalib*1.2);
        funcCalib.SetRange(peakCalib*0.7,peakCalib*1.3);
        funcCalib.SetParameter(1,peakCalib);
        hist->Fit(&funcCalib, "0QR");
        calibValue[i]=funcCalib.GetParameter(1);
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
      if(!iter->GetRecoPulseFlags().none())
      {
        if(!iter->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::noCalibConstPulseArea) &&
           !iter->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::noCalibConstPulseHeight)) continue;
      }

      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      int channelIndex=barIndex*CRVId::nChanPerBar+SiPM;
      _calibHistsPulseArea.at(channelIndex)->Fill(iter->GetPulseBeta()*iter->GetPulseHeight()*TMath::E());
      _calibHistsPulseHeight.at(channelIndex)->Fill(iter->GetPulseHeight());
    }
  }

} // end namespace mu2e

using mu2e::CrvCalibration;
DEFINE_ART_MODULE(CrvCalibration)
