//
// A module to convert simulated CRV (Wideband) hits into a ROOT file
// with a structure similar to the one used for the Wideband standalone analysis.
//
//
// Original Author: Ralf Ehrlich

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CRVConditions/inc/CRVCalib.hh"

#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"

#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"

#include <string>

#include <TMath.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFitResult.h>

namespace mu2e
{
  class CrvTimingStudies : public art::EDAnalyzer
  {
    public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config
    {
      fhicl::Atom<std::string> crvRecoPulsesModuleLabel{ Name("crvRecoPulseModuleLabel"), Comment("CrvRecoPulse Label")};
      fhicl::Atom<double>      PEthreshold{ Name("PEthreshold"), Comment("PE threshold")};
      fhicl::Atom<bool>        removeTimeOffsets{ Name("removeTimeOffsets"), Comment("remove time offsets added by reco")};
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CrvTimingStudies(const Parameters& conf);
    ~CrvTimingStudies() override;
    void analyze(const art::Event& e) override;
    void beginJob() override;
    void endJob() override;

    private:
    std::string   _crvRecoPulsesModuleLabel;
    double        _PEthreshold;
    bool          _removeTimeOffsets;

    std::map<std::pair<int,int>,TH1F*>  _histTimeDiffs;

    ProditionsHandle<CRVOrdinal> _crvChannelMap_h;
    ProditionsHandle<CRVCalib>   _calib_h;
  };

  CrvTimingStudies::CrvTimingStudies(const Parameters& conf) :
    art::EDAnalyzer{conf},
    _crvRecoPulsesModuleLabel(conf().crvRecoPulsesModuleLabel()),
    _PEthreshold(conf().PEthreshold()),
    _removeTimeOffsets(conf().removeTimeOffsets())
  {
  }

  CrvTimingStudies::~CrvTimingStudies()
  {
  }

  void CrvTimingStudies::beginJob()
  {
  }

  void CrvTimingStudies::endJob()
  {
  }

  void CrvTimingStudies::analyze(const art::Event& event)
  {
    art::Handle<CrvRecoPulseCollection> crvRecoPulseCollection;
    if(!event.getByLabel(_crvRecoPulsesModuleLabel,"",crvRecoPulseCollection)) return;

    auto const& crvChannelMap = _crvChannelMap_h.get(event.id());
    auto const& calib = _calib_h.get(event.id());

    static bool firstEvent=true;
    if(firstEvent)
    {
      firstEvent=false;

      //store channel map, pedestals and calibration constants in the file,
      //so that it can later be used to write a full calibration set
      //of pedestals, calib consts, and time offsets

      GeomHandle<CosmicRayShield> CRS;
      const std::vector<std::shared_ptr<CRSScintillatorBar> > &counters = CRS->getAllCRSScintillatorBars();
      art::ServiceHandle<art::TFileService> tfs;

      TTree *treeChannelMap = tfs->make<TTree>("channelMap","channelMap");
      size_t channel;
      int    roc, feb, febChannel;
      treeChannelMap->Branch("channel", &channel);
      treeChannelMap->Branch("roc", &roc);
      treeChannelMap->Branch("feb", &feb);
      treeChannelMap->Branch("febChannel", &febChannel);

      for(channel=0; channel<counters.size()*CRVId::nChanPerBar; ++channel)
      {
        if(!crvChannelMap.onlineExists(channel)) continue;
        CRVROC onlineChannel  = crvChannelMap.online(channel);
        roc         = onlineChannel.ROC();
        feb         = onlineChannel.FEB();
        febChannel  = onlineChannel.FEBchannel();
        treeChannelMap->Fill();
      }

      TTree *treeCalib = tfs->make<TTree>("crvCalib","crvCalib");
      double pedestal, calibPulseHeight, calibPulseArea;
      treeCalib->Branch("channel", &channel);
      treeCalib->Branch("pedestal", &pedestal);
      treeCalib->Branch("calibPulseHeight", &calibPulseHeight);
      treeCalib->Branch("calibPulseArea", &calibPulseArea);

      for(channel=0; channel<counters.size()*CRVId::nChanPerBar; ++channel)
      {
        pedestal = calib.pedestal(channel);
        calibPulseHeight = calib.pulseHeight(channel);
        calibPulseArea = calib.pulseArea(channel);
        treeCalib->Fill();
      }
    }

    int previousOfflineChannel=-1;
    std::map<uint16_t,std::vector<double> > fpgaTimes;
    for(auto crvRecoPulse=crvRecoPulseCollection->begin(); crvRecoPulse!=crvRecoPulseCollection->end(); ++crvRecoPulse)
    {
      if(crvRecoPulse->GetPEs()<_PEthreshold) continue;  //ignore pulses below a threshold
      if(crvRecoPulse->GetRecoPulseFlags().test(CrvRecoPulseFlagEnums::failedFit)) continue;  //ignore pulses with a failed fit

      auto barIndex         = crvRecoPulse->GetScintillatorBarIndex();
      int  SiPM             = crvRecoPulse->GetSiPMNumber();
      double recoPulseTime  = crvRecoPulse->GetPulseTime();

      uint16_t offlineChannel = barIndex.asUint()*CRVId::nChanPerBar + SiPM;
      if(offlineChannel==previousOfflineChannel) continue;
      previousOfflineChannel=offlineChannel;
      CRVROC   onlineChannel  = crvChannelMap.online(offlineChannel);
      uint16_t ROC            = onlineChannel.ROC();
      uint16_t feb            = onlineChannel.FEB();
      uint16_t febChannel     = onlineChannel.FEBchannel();

      uint16_t fpgaIndex      = ((ROC-1)*CRVId::nFEBPerROC*CRVId::nChanPerFEB+(feb-1)*CRVId::nChanPerFEB+febChannel)/(CRVId::nChanPerFEB/CRVId::nFPGAPerFEB);

      if(_removeTimeOffsets) //remove time offsets introduced during reconstruction to see the "pure" time differences, e.g. to generate new time calibration tables
      {
        double timeOffset = calib.timeOffset(offlineChannel);
        recoPulseTime-=timeOffset;
      }

      fpgaTimes[fpgaIndex].push_back(recoPulseTime);
    }

    std::map<uint16_t,double> fpgaAverageTimes;
    for(auto fpga=fpgaTimes.begin(); fpga!=fpgaTimes.end(); ++fpga)
    {
      const auto &times=fpga->second;
      double averageTime=0;
      for(size_t i=0; i<times.size(); ++i) averageTime+=times.at(i);
      averageTime/=times.size();
      fpgaAverageTimes[fpga->first]=averageTime;
    }

    //compare time differences between different FPGAs
    for(auto fpga1=fpgaAverageTimes.begin(); fpga1!=fpgaAverageTimes.end(); fpga1++)
    for(auto fpga2=fpgaAverageTimes.begin(); fpga2!=fpgaAverageTimes.end(); fpga2++)
    {
      if(fpga1->first>=fpga2->first) continue; //don't compare with itself (=) and avoid comparing the same FPGAs twice (>)

      std::pair<int,int> histIndex(fpga1->first,fpga2->first);
      if(_histTimeDiffs.find(histIndex)==_histTimeDiffs.end())
      {
        art::ServiceHandle<art::TFileService> tfs;
        _histTimeDiffs[histIndex] = tfs->make<TH1F>(Form("fpgaTimeDiff_%i_%i",fpga1->first,fpga2->first),
                                                    Form("Time Diffs between FGPAs %i and %i;time difference [ns];Counts",fpga1->first,fpga2->first),
                                                    300,-150,150);
      }

      double timeDiff=fpga1->second-fpga2->second;
      _histTimeDiffs[histIndex]->Fill(timeDiff);
    }

  } // end analyze

} // end namespace mu2e

using mu2e::CrvTimingStudies;
DEFINE_ART_MODULE(CrvTimingStudies)
