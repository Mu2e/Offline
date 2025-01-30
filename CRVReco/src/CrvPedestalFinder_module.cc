//
// A module that finds the pedestals for Crv channels
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/DataProducts/inc/CRVId.hh"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include <TH1F.h>
#include <TF1.h>

namespace mu2e
{
  class CrvPedestalFinder : public art::EDAnalyzer
  {

    public:
    struct Config
    {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> crvDigiModuleLabel{Name("crvDigiModuleLabel"), Comment("module label for CrvDigis")};
      fhicl::Atom<double>      maxADCspread{Name("maxADCspread"), Comment("maximum spread of ADC values within a waveform to be considered for the pedestal")};
      fhicl::Atom<std::string> tmpDBfileName{Name("tmpDBfileName"), Comment("name of the tmp. DB file name for the pedestals")};
    };

    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit CrvPedestalFinder(const Parameters& config);
    void analyze(const art::Event& e);
    void beginRun(const art::Run&);
    void endRun(const art::Run&);

    private:
    std::string        _crvDigiModuleLabel;
    double             _maxADCspread;
    std::string        _tmpDBfileName;
    std::vector<TH1F*> _pedestalHists;

    ProditionsHandle<CRVCalib> _calib_h;

    std::vector<double> _timeOffsets;
  };


  CrvPedestalFinder::CrvPedestalFinder(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _crvDigiModuleLabel(conf().crvDigiModuleLabel()), _maxADCspread(conf().maxADCspread()), _tmpDBfileName(conf().tmpDBfileName())
  {
  }

  void CrvPedestalFinder::beginRun(const art::Run&)
  {
    GeomHandle<CosmicRayShield> CRS;
    const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
    _pedestalHists.reserve(counters.size()*CRVId::nChanPerBar);
    _timeOffsets.resize(counters.size()*CRVId::nChanPerBar);

    art::ServiceHandle<art::TFileService> tfs;
    for(size_t barIndex=0; barIndex<counters.size(); ++barIndex)
    {
      for(size_t SiPM=0; SiPM<CRVId::nChanPerBar; ++SiPM)
      {
        //produce histograms also for non-existing channels to get a continuously running index
        size_t channelIndex=barIndex*4 + SiPM;
        _pedestalHists.emplace_back(tfs->make<TH1F>(Form("crvPedestalHist_%lu",channelIndex),
                                                    Form("crvPedestalHist_%lu",channelIndex),
                                                    200,-50,150));   //TODO: needs to be only between -50 and +50, but Offline currently sets the pedestal at +100
        _timeOffsets[channelIndex]=0;
      }
    }
  }

  void CrvPedestalFinder::endRun(const art::Run&)
  {
    TF1 funcPedestal("f0", "gaus");

    std::ofstream outputFile;
    outputFile.open(_tmpDBfileName);
    outputFile<<"TABLE CRVSiPM"<<std::endl;
    outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

    for(size_t channel=0; channel<_pedestalHists.size(); ++channel)
    {
      TH1F *hist=_pedestalHists.at(channel);
      if(hist->GetEntries()<100) //not enough data
      {
        outputFile<<channel<<","<<0<<",-1,-1"<<std::endl;
        continue;
      }

      int n=hist->GetNbinsX();
      double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
      if(overflow/((double)hist->GetEntries())>0.1) //too much underflow/overflow. something may be wrong.
      {
        outputFile<<channel<<","<<0<<",-1,-1"<<std::endl;
        continue;
      }

      int maxbinPedestal = hist->GetMaximumBin();
      double peakPedestal = hist->GetBinCenter(maxbinPedestal);
      funcPedestal.SetRange(peakPedestal-4,peakPedestal+4);
      funcPedestal.SetParameter(1,peakPedestal);
      hist->Fit(&funcPedestal, "0QR");
      outputFile<<channel<<","<<funcPedestal.GetParameter(1)<<",-1,-1"<<std::endl;  //only print out pedestal values
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

  void CrvPedestalFinder::analyze(const art::Event& event)
  {
    art::Handle<CrvDigiCollection> crvDigiCollection;
    if(!event.getByLabel(_crvDigiModuleLabel,"NZS",crvDigiCollection)) return;

    //find time offsets from first event
    //need to assume that this is only used for calibration runs where the time offsets stay constant over the entire run
    static bool first=true;
    if(first)
    {
      first=false;
      auto const& calib = _calib_h.get(event.id());
      for(size_t channelIndex=0; channelIndex<_timeOffsets.size(); ++channelIndex)
      {
        _timeOffsets[channelIndex] = calib.timeOffset(channelIndex);
      }
    }

    for(auto iter=crvDigiCollection->begin(); iter!=crvDigiCollection->end(); ++iter)
    {
      auto minmaxTest = std::minmax_element(iter->GetADCs().begin(),iter->GetADCs().end());
      if(*minmaxTest.second-*minmaxTest.first>_maxADCspread) continue;  //ignore waveforms that fluctuate too much.
                                                                        //they probably contain a signal or noise pulse

      int barIndex = iter->GetScintillatorBarIndex().asInt();
      int SiPM = iter->GetSiPMNumber();
      int channelIndex=barIndex*CRVId::nChanPerBar+SiPM;
      auto hist = _pedestalHists.at(channelIndex);

      for(size_t i=0; i<iter->GetADCs().size(); ++i)
      {
        hist->Fill(iter->GetADCs()[i]);
      }
    }
  }

} // end namespace mu2e

using mu2e::CrvPedestalFinder;
DEFINE_ART_MODULE(CrvPedestalFinder)
