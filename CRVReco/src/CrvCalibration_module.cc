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
      fhicl::Atom<std::string> tmpDBfileName{Name("tmpDBfileName"), Comment("name of the tmp. DB file name for the pedestals")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvCalibration(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);
    void endRun(const art::Run&);

    private:
    std::string        _crvRecoPulsesModuleLabel;
    std::string        _tmpDBfileName;
    std::vector<TH1F*> _calibHistsPulseArea;
    std::vector<TH1F*> _calibHistsPulseHeight;

    ProditionsHandle<CRVCalib> _calib_h;

    std::vector<double> _pedestals;
    std::vector<double> _timeOffsets;
  };


  CrvCalibration::CrvCalibration(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()), _tmpDBfileName(conf().tmpDBfileName())
  {
  }

  void CrvCalibration::beginRun(const art::Run&)
  {
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
        size_t channelIndex=barIndex*4 + SiPM;
        _calibHistsPulseArea.emplace_back(tfs->make<TH1F>(Form("crvCalibrationHistPulseArea_%lu",channelIndex),
                                                 Form("crvCalibrationHistPulseArea_%lu",channelIndex),
                                                 150,0,3000));
        _calibHistsPulseHeight.emplace_back(tfs->make<TH1F>(Form("crvCalibrationHistPulseHeight_%lu",channelIndex),
                                                 Form("crvCalibrationHistPulseHeight_%lu",channelIndex),
                                                 150,0,75));
        _pedestals[channelIndex]=0;
        _timeOffsets[channelIndex]=0;
      }
    }
  }

  void CrvCalibration::endRun(const art::Run&)
  {
    TF1 funcCalib("f0", "gaus");

    std::ofstream outputFile;
    outputFile.open(_tmpDBfileName);
    outputFile<<"TABLE CRVSiPM"<<std::endl;
    outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

    for(size_t channel=0; channel<_pedestals.size(); ++channel)
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

        int n=hist->GetNbinsX();
        double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
        if(overflow/((double)hist->GetEntries())>0.1) //too much underflow/overflow. something may be wrong.
        {
          calibValue[i]=-1;
          continue;
        }

        int maxbinCalib = hist->GetMaximumBin();
        double peakCalib = hist->GetBinCenter(maxbinCalib);
        funcCalib.SetRange(peakCalib-4,peakCalib+4);
        funcCalib.SetParameter(1,peakCalib);
        hist->Fit(&funcCalib, "NQR");
        calibValue[i]=funcCalib.GetParameter(1);
      }

      outputFile<<channel<<","<<_pedestals.at(channel)<<","<<calibValue[0]<<","<<calibValue[1]<<std::endl;
    }

    outputFile<<std::endl;

    //CRVTime table
    outputFile<<"TABLE CRVTime"<<std::endl;
    outputFile<<"#channel, timeOffset"<<std::endl;
    for(size_t channel=0; channel<_timeOffsets.size(); ++channel)
    {
      outputFile<<channel<<","<<_timeOffsets.at(channel)<<std::endl;
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
    }

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
